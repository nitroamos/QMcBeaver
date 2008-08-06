#!/usr/bin/perl
use POSIX;

#alphabet first, numerical second
sub a1n2 {
    my @adata = split/&/,$a;
    my @bdata = split/&/,$b;
    $bdata[1] <=> $adata[1];
    if($adata[0] eq $bdata[0]){
	if($adata[3] eq $bdata[3]) {
	    $bdata[1] cmp $adata[1];
	} else {
	    $adata[3] <=> $bdata[3];
	}
    } else {
	$bdata[0] cmp $adata[0];
    }
}

sub byenergy {
    my @adata = split/&/,$a;
    my @bdata = split/&/,$b;
    if($adata[0] != $bdata[0]){
	$bdata[0] <=> $adata[0];
    } else {
	$bdata[1] <=> $adata[1];
    }
}

sub bydt {
    my @adata = split/&/,$a;
    my @bdata = split/&/,$b;
    if($adata[1] != $bdata[1]){
	#compare dt
	$bdata[1] <=> $adata[1];
    } else {
	#compare energies
	$bdata[0] <=> $adata[0];
    }
}

sub gcf {
    my ($x, $y) = @_;
    ($x, $y) = ($y, $x % $y) while $y;
    return $x;
}

sub getEnergyWError {
    my ($nrg, $err) = @_;
    my $str = "";
    if(abs($err) == 0){	
	$str = "$nrg";
    } else {
	my $d = 1-int(floor(log($err)/log(10.0)));
	my $energy = floor( $nrg * pow(10.0,$d) + 0.5) / pow(10.0,$d);
	$str = sprintf "%.*f",$d,$energy;
	my $error  = floor( $err * pow(10.0,$d) + 0.5);
	$str = "$str($error)";
    }
    #printf("nrg=%10.5f err=%10.5f d=%3i energy=%20f str=%s\n",$nrg,$err,$d,$energy,$str);
    return $str;
}
sub getFormula {
    my ($one, $two) = @_;
    my @od = split/&/,$one;
    my @td = split/&/,$two;

    return (0,0) if($od[0] == $td[0] || #compare energies
		    $od[1] != $td[1] || #compare dt
		    $od[4] != $td[4] || #compare num walkers
		    $od[3] ne $td[3]);  #compare jastrows

    $or = $od[2];
    $tr = $td[2];
    $factor = 100;

    while($factor != 1){
	$factor = gcf($or,$tr);
	#print "gcf($or,$tr) = $factor\n";
	$or /= $factor;
	$tr /= $factor;
    }

    #print "ratio = $ratio : intr = $intr : intdiff = $intdiff\n";
    if($or == int($or) && $tr == int($tr) &&
       $or < 10 && $tr < 10)
    {
	return ($tr, $or);
    } else {
	return (0, 0);
    }
}

sub estimateTimeToFinish
{
    my ($outfile) = @_;
    return 0 if(!(-e $outfile));

    @newsteps = `grep "new steps" $outfile`;

    my $totalSteps = 0;
    if($#newsteps  < 0){
	my $infile = substr($outfile,0,-4);
	$infile .= ".ckmf";
	open(CKMFFILE,"$infile");
	while(<CKMFFILE>){
	    if($_ =~ m/^\s*max_time_steps\s*$/){
		$_ = <CKMFFILE>;
		chomp;
		my @line = split/[ ]+/;
		$totalSteps = $line[1];
		last;
	    }
	}
    } else {
	$totalSteps = (split/\s+/,$newsteps[$#newsteps])[12];
    }


    @itertime = `grep "Average iterations per hour:" $outfile`;
    return 0 if($#itertime < 0);
    my $shift = $#itertime;
    #the correlated sampling phase runs faster per iteration, and we assume that we're currently
    #in the longer phase, so we want to look back 2 iterations
    $shift -= 1 if($shift > 0);
    my $itersPerHour = (split/\s+/,$itertime[$shift])[4];

    my $curIter = (split/\s+/,`tail -n 1 $outfile`)[1];

    my $est = ($totalSteps - $curIter) / $itersPerHour;

    my $hrs = int($est);
    my $mns = int(($est - $hrs)*60.0 + 0.5);
    if($mns < 10){
	$mns = "0$mns";
    }
    #print "$est = hrs = $hrs mns = $mns\n";
    #print "totalSteps = $totalSteps curiter = $curIter itertime = $itersPerHour est = $est\n";
    #return $est;
    return "${hrs}:${mns}";
}

sub getCKMFHeader
{
    my $header = sprintf "%-30s %2s %3s %11s %-7s %-7s %1s%1s%1s%1s %-15s %8s %8s\n",
    "Name","  ","NW","EQ/Steps","dt","nci:nbf","O","L","C","3","HF","Age","Size";
    return $header;
}

sub getCKMFSummary
{
    my ($ckmf) = @_;

    $base = substr($ckmf,0,-5);
    my $dirname   = `dirname $base`;
    chomp($dirname);
    my $shortbase = `basename $base`;

    if($dirname eq "."){

    } else {
	my $nextbase = `basename $dirname`;
	chomp($nextbase);
	$shortbase = "$nextbase/$shortbase";
    }

    chomp ($shortbase);
    my $curTime = qx! date +%s !; 
    open (CKMFFILE, "$ckmf");

    while(<CKMFFILE> !~ /&flags/){};
    $rt      = "";
    $numbf   = 0;
    $numci   = 0;
    $hfe     = "";
    $nw      = 0;
    $dt      = 0;
    $opt     = 0;
    $optl    = 0;
    $optci   = 0;
    $opt3    = 0;
    $steps   = 0;
    $eqsteps = 0;
    $iseed   = 0;

    while(<CKMFFILE>){
	if($_ =~ m/^\s*run_type\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $rt = $line[1];
	}
	if($_ =~ m/^\s*energy\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $hfe = $line[1];
	}
	if($_ =~ m/^\s*nbasisfunc\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $numbf = $line[1];
	}
	if($_ =~ m/^\s*ndeterminants\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $numci = $line[1];
	}
	if($_ =~ m/^\s*dt\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $dt = $line[1];
	}
	if($_ =~ m/^\s*number_of_walkers\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $nw = $line[1];
	}
	if($_ =~ m/^\s*optimize_Psi\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $opt = $line[1];
	}
	if($_ =~ m/^\s*optimize_L\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $optl = $line[1];
	}
	if($_ =~ m/^\s*optimize_CI\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $optci = $line[1];
	}
	if($_ =~ m/^\s*optimize_NEE_Jastrows\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $opt3 = $line[1];
	}
	if($_ =~ m/^\s*max_time_steps\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $steps = $line[1];
	}
	if($_ =~ m/^\s*equilibration_steps\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $eqsteps = $line[1];
	}
	if($_ =~ m/^\s*iseed\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $iseed = $line[1];
	}

	#any other interesting parameters?
	#if($ARGV[0] != "" && $_ =~ m/$ARGV[0]/ && $_ !~ /\#/){
	#    $name = $_;
	#    chomp $name;
	#    $_ = <CKMFFILE>;
	#    chomp;
	#    my @line = split/[ ]+/;
	#    $val = $line[1];
	#    $searchdata .= sprintf "%20s: %30s = %30s\n", "", $name, $val;
	#}
	if($_ =~ m/&geometry$/){
	    last;
	}
    }

    if($rt eq "variational"){
	$rt = "v";
    } elsif($rt = "diffusion"){
	$rt = "d";
    }

    my $outModTime = "";
    my $outSize    = "";
    my $ovData     = "";
    my $failed     = 0;

    if(-e "$base.out"){
	my $data = `/bin/ls -lh --time-style=+%s $base.out`;
	my @list = split/ +/,$data;

	$outSize = $list[4];
	$outModTime = $curTime - $list[5];
	$char = " ";

	if($outModTime > 3600){
	    $outModTime /= 3600;
	    $char = "h";
	    if($outModTime > 24){
		$outModTime /= 24;
		$char = "d";
	    }
	}
	if($char eq " "){
	    $outModTime = sprintf "%5.0f $char", $outModTime;
	} else {
	    $outModTime = sprintf "%5.1f $char", $outModTime;
	}
	#$outModTime .= sprintf " %3s", $list[5];
	#$outModTime .= sprintf " %2s", $list[6];
	#$outModTime .= sprintf " %5s", $list[7];
    }

    if(-e "$base.out" && $opt){
	$data = `grep failed $base.out`;
	$failed = 1 if(length($data) > 0);

	if($failed == 1){
	    $outModTime .= "*";
	} else {
	    $outModTime .= " ";
	}

	@list = `grep "Objective Value" $base.out`;
	$ovData = $list[$#list];
	chomp($ovData);

	@list = split/[ =]+/,$ovData;
	$ovData = "";
	if($#list == 9){
	    $ovData .= sprintf "%2i",$list[1];
	    $ovData .= sprintf " %15.10f",$list[5];
	    $ovData .= sprintf " %8.5f",$list[7];
	    $ovData .= sprintf " %10s",$list[9];
	}

	@newsteps = `grep "new steps" $base.out`;
	my $curSteps = $steps;
	if($#newsteps  >= 0){	   
	    $curSteps = (split/\s+/,$newsteps[$#newsteps])[12];	
	}
	$ovData .= sprintf " %10s",$curSteps;
    }

    $steps_str = "";
    if($steps >= 1000*1000*1000){
	$steps /= int(1000*1000*1000);
	$steps_str = "${steps}B";
    } elsif($steps >= 1000*1000){
	$steps /= int(1000*1000);
	$steps_str = "${steps}M";
    } elsif($steps >= 1000){
	$steps /= int(1000);
	$steps_str = "${steps}K";
    } else {
	$steps_str = "$steps";
    }

    $eqsteps_str = "";
    if($eqsteps >= 1000*1000*1000){
	$eqsteps /= int(1000*1000*1000);
	$eqsteps_str = "${eqsteps}B";
    } elsif($eqsteps >= 1000*1000){
	$eqsteps /= int(1000*1000);
	$eqsteps_str = "${eqsteps}M";
    } elsif($eqsteps >= 1000){
	$eqsteps /= int(1000);
	$eqsteps_str = "${eqsteps}K";
    } else {
	$eqsteps_str = "$eqsteps";
    }

    my $oneliner = "";
    $oneliner .= sprintf "%-30s %2s %3i %5s/%-5s %-7s %3i:%-3s %1i%1i%1i%1i %-15s",
    $shortbase,$rt,
    $nw,$eqsteps_str,$steps_str,
    $dt,
    ${numci},${numbf},$opt,$optl,$optci,$opt3,$hfe;
    $oneliner .= sprintf " %10s", $outModTime;
    $oneliner .= sprintf " %7s", $outSize;
    $oneliner .= sprintf " %50s", $ovData;
    if($iseed != 0){
	$oneliner .= sprintf " iseed = $iseed";
    }
    $oneliner .= sprintf "\n";

    close(CKMFFILE);
    return $oneliner;
}

sub getEnergies
{
    my ($filename, $energies) = @_;
    open(FILE,"$filename");

    $more = 1;
    while(<FILE>){
	$sampleclock = (split/[ ]+/)[8] if(/Average microseconds per sample per num initial walkers/);
	$sampleVar = (split/[ ]+/)[3] if(/Sample variance/ && $sampleVar == 0);	    

	#this is to avoid processing warnings
	next if( $_ =~ /[=:]/ && $_ !~ /Results/);

	chomp;
	@data = split/[ ]+/;

        #this is the number of data elements per line
	#it can have the letter 'e' or 'E' since scientific notation uses them
	if($#data >= 8 && $_ !~ /[A-DF-Za-df-z]+/ && $more){
	    $counter++;
	    $iteration   = $data[1];
	    $eavg        = $data[2];
	    $estd        = $data[3];
	    if(abs($eavg) > 1e-10){
		push(@$energies,$eavg);
	    }
	} elsif(/Results/) {
	    #$more = 0;
	}
    }    
    close(FILE);
}

# this function will fill in the files array with
# all files that have the extension ext. there are
# a few known directories it will not descend into
sub getFileList
{
    my ($ext, $files) = @_;

    #this will scan through all the subdirectories in the $files array looking for $ext files
    my $clean = 0;
    my $loops = 0;
    while($clean == 0)
    {
	$loops++;
	$clean = 1;
	my @newfiles;
	
	for(my $index=0; $index<=$#$files; $index++)
	{
	    my $cur = ${@$files}[$index];
	    chomp($cur);

	    #there are some obvious directories we don't need to search.
	    #we also don't look in folders that end in 'hide', unless it was specified on the command line
	    if(-d $cur && $cur !~ /src$/ && $cur !~ /bin$/ && $cur !~ /include$/ && 
	       ($cur !~ /hide$/ || $loops <= 1))
	    {
		my @list = `ls $cur`;
		foreach $item (@list)
		{
		    #we have a directory in the list, so we're going to need to loop again
		    $clean = 0;
		    chomp($item);
		    if($cur eq ".")
		    {
			push(@newfiles,"$item");
		    } else {
			push(@newfiles,"$cur/$item");
		    }
		}	    
	    } elsif($cur =~ /$ext$/ && $cur !~ /.step[\d]+./ && $cur !~ /.opt[\d]+./){
		push(@newfiles,$cur);
	    }
	}
	@$files = @newfiles;

	if($loops > 8)
	{
	    print "Stopping recursion at $loops.\n";
	}
    }
}
1;
