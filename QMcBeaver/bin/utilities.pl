#!/usr/bin/perl
use POSIX;

sub areComparable {
# This function is used by the code to see if two calculations can be compared.
# The script will generate output comparing each result against all other results,
# which add up to quite a few comparisons, most of which are actually meaningless.
# So if they're meaningless, then return 0. You'll probably want to edit this function
# to choose your own comparisons.
#
# The input is from summary.pl, where each a key is created for each calculation:
# my $key = "$refE&$dt&$numbf&$numjw&$nw&$numci&$numor&$oepi&$short"; 
#
    my ($one, $two) = @_;
    my @od = split/&/,$one;
    my @td = split/&/,$two;

    return 0 if($od[0] == $td[0] || #compare energies
		$od[1] != $td[1] || #compare dt
		$od[4] != $td[4] || #compare num walkers
		$od[7] != $td[7]); #compare oepi

    #make sure the jastrows are comparable
    return 0 if($od[3] =~ /44/ && $td[3] !~ /44/);
    return 0 if($od[3] !~ /44/ && $td[3] =~ /44/);
    
    #the files are named something like awt0p2, so extract the letter after the 0 (or 4),
    #p in this case, and make sure they match
    my $oType = "";
    my $tType = "";
    $oType = $1 if($od[8] =~ /t\d(\w)/);
    $tType = $1 if($td[8] =~ /t\d(\w)/);
    #this probably needs to be turned off for atomization energies
    return 0 if($oType ne $tType);

    #make sure the last number in the file matches
    #This only makes a difference if we didn't average over the results.
    my $oLast = "";
    my $tLast = "";
    $oLast = $1 if($od[8] =~ /([\d\.]+)$/);
    $tLast = $1 if($td[8] =~ /([\d\.]+)$/);
    #return 0 if($oLast ne $tLast);

    return 1;
}

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

sub a2n3 {
    my @adata = split/&/,$a;
    my @bdata = split/&/,$b;
    if($adata[0] eq $bdata[0]){
	if($adata[5] eq $bdata[5]){
	    #sort by opt iter
	    $bdata[1] <=> $adata[1];
	} else {
	    #sort by reference energy
	    $adata[5] cmp $bdata[5];
	}
    } else {
	#Sort by jastrow type (e.g. s, t, UC, etc)
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
    my ($a, $b, $c, $orbFilter) = @_;
    my $am = $c;
    my $bm = $b;
    my $cm = $a;

    my $factor = 100;
    while($factor != 1){
	$factor = gcf($am,$cm);
	#print "gcf($ar,$cr) = $factor\n";
	$am /= $factor;
	$cm /= $factor;
    }

    if($am == int($am) && $cm == int($cm) &&
       $am < 10 && $cm < 10)
    {
	#return (0,0) if($ar*$cd[6] != $cr*$ad[6] &&
	#		$arbFilter);
	#print "($a, 0, $c) => ($am, 0, $cm)\n";
	return ($am, 0, $cm);
    }

    my $maxF = 6;
    for($am=1; $am <= $maxF; $am+=1){
	for($bm=1; $bm <= $maxF; $bm+=1){
	    for($cm=1; $cm <= $maxF; $cm+=1){
		if( $am*$a + $bm*$b == $cm*$c){
		    #print "($a, $b, $c) => ($am, $bm, $cm)\n";
		    return ($am,$bm,$cm);
		}
	    }
	}
    }    
    
    return (0, 0, 0);
}

sub getFileAge
{
    my ($file,$abstime) = @_;
    my $curTime = qx! date +%s !; 
    my $data = `/bin/ls -lh --time-style=+%s $file`;
    my @list = split/ +/,$data;
    $outSize = $list[4];
    my $outModTime = $curTime - $list[5];

    return $outModTime if($abstime == 1);
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
    return $outModTime;
}

sub estimateTimeToFinish
{
    my ($outfile, $time) = @_;
    return 0 if(!(-e $outfile));
    my $base = substr($outfile,0,-4);
    @newsteps = `grep "new steps" $outfile`;

    my $equilSteps = 0;
    my $totalSteps = 0;
    if($#newsteps  < 0){
	open(CKMFFILE,"${base}.ckmf");
	while(<CKMFFILE>){
	    if($_ =~ m/^\s*max_time_steps\s*$/){
		$_ = <CKMFFILE>;
		chomp;
		my @line = split/[ ]+/;
		$totalSteps += $line[1];
	    }
	    if($_ =~ m/^\s*equilibration_steps\s*$/){
		$_ = <CKMFFILE>;
		chomp;
		my @line = split/[ ]+/;
		$equilSteps = $line[1];
		#$totalSteps += $line[1];
	    }
	}
    } else {
	$totalSteps = (split/\s+/,$newsteps[$#newsteps])[12];
    }

    @itertime = `grep "Average iterations per hour:" $outfile`;
    my $curIter = (split/\s+/,`tail -n 1 ${base}.qmc`)[1];
    $curIter += $equilSteps if($curIter <= 0);

    my $itersPerHour = 0;
    if($#itertime < 0 && $time != 0){
	$itersPerHour = $curIter / $time;
	$itersPerHour *= 3600;
    } elsif($#itertime >= 0) {
	my $shift = $#itertime;
	#the correlated sampling phase runs faster per iteration, and we assume that we're currently
	#in the longer phase, so we want to look back 2 iterations
	$shift -= 1 if($shift > 0);
	$itersPerHour = (split/\s+/,$itertime[$shift])[4];       
    } else {
	return "0:0";
    }
    return "0" if($itersPerHour == 0);
    my $est = ($totalSteps - $curIter) / $itersPerHour;
    my $hrs = int($est);
    my $mns = int(($est - $hrs)*60.0 + 0.5);
    if($mns < 10){
	$mns = "0$mns";
    }
    #print "$est => hrs = $hrs mns = $mns\n";
    #print "totalSteps = $totalSteps curiter = $curIter itertime = $itersPerHour est = $est\n";
    return "${hrs}:${mns}";
}

sub getOPTHeader
{
    return "IDUE3L";
}
sub getCKMFHeader
{
    my $header = sprintf "%69s%5s\n",""," CUUN";

    $header .= sprintf "%-30s %2s O %3s %11s %1s %-7s %-7s %6s %-15s %8s %8s\n",
    "Name","  ","NW","EQ/Steps","e","dt","nci:nbf",
    getOPTHeader(),
    "HF","Age","Size";
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
    open (CKMFFILE, "$ckmf");

    while(<CKMFFILE> !~ /&flags/){};
    $rt      = "";
    $numbf   = 0;
    $numci   = 0;
    $hfe     = "";
    $nw      = 0;
    $dt      = 0;
    $steps   = 0;
    $eqsteps = 0;
    $iseed   = 0;
    $oepi    = 0;

    $opt     = 0;
    $optl    = 0;
    $optci   = 0;
    $opt3    = 0;
    $optUD   = 0;
    $optUU   = 0;
    $optNE   = 0;
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
	if($_ =~ m/^\s*one_e_per_iter\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $oepi = $line[1];
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
	if($_ =~ m/^\s*optimize_EE_Jastrows\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $optUU = $line[1];
	    $optUD = $line[1];
	}
	if($_ =~ m/^\s*optimize_EN_Jastrows\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $optNE = $line[1];
	}
	if($_ =~ m/^\s*optimize_UD_Jastrows\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $optUD = $line[1];
	}
	if($_ =~ m/^\s*optimize_UU_Jastrows\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $optUU = $line[1];
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
	$outModTime = getFileAge("$base.out",0);
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
	$steps_str = sprintf "%2.1fB",${steps}; 
    } elsif($steps >= 1000*1000){
	$steps /= int(1000*1000);
	$steps_str = sprintf "%2.1fM",${steps};
    } elsif($steps >= 1000){
	$steps /= int(1000);
	$steps_str = sprintf "%2.1fK",${steps}; 
    } else {
	$steps_str = "$steps";
    }

    $eqsteps_str = "";
    if($eqsteps >= 1000*1000*1000){
	$eqsteps /= int(1000*1000*1000);
	$eqsteps_str = sprintf "%2.1fB",${eqsteps};
    } elsif($eqsteps >= 1000*1000){
	$eqsteps /= int(1000*1000);
	$eqsteps_str = sprintf "%2.1fM",${eqsteps};  
    } elsif($eqsteps >= 1000){
	$eqsteps /= int(1000);
	$eqsteps_str = sprintf "%2.1fK",${eqsteps};  
    } else {
	$eqsteps_str = "$eqsteps";
    }

    my $oneliner = "";
    $oneliner .= sprintf "%-30s %2s %1i %3i %5s/%-5s %1s %-7s %3i:%-3s %1i%1i%1i%1i%1i%1i %-15s",
    $shortbase,$rt,$opt,
    $nw,$eqsteps_str,$steps_str,
    $oepi,$dt,
    ${numci},${numbf},
    $optci,$optUD,$optUU,$optNE,$opt3,$optl,$hfe;
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
		#turn all // in file paths to just one /
		$cur =~ s/\/\//\//;
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
