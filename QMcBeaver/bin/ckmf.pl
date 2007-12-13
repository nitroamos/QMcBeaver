#!/usr/bin/perl

my @files;
if($#files < 0)
{
    push(@files,".");
    #if the first argument is a plotfile.dat, then we'll add new data to it
    #otherwise, we'll create a new plotfile.dat
    
    #print "Usage: qmc_convergence_graph.pl {plotfile.dat || output1.out} [output2.out ...]\n";
    #print "Will produce a graph of energy vs num samples or iterations.\n";
    #die;
}

#this will scan through all the subdirectories looking for .out files
my $clean = 0;
my $loops = 0;
while($clean == 0)
{
    $loops++;
    $clean = 1;
    my @newfiles;

    for(my $index=0; $index<=$#files; $index++)
    {
	my $cur = $files[$index];
	chomp($cur);
	if(-d $cur && $cur !~ /src$/ && $cur !~ /bin$/ && $cur !~ /include$/ && $cur !~ /hide$/)
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
	} elsif($cur =~ /.ckmf$/ && $cur !~ /.step[\d]+./ && $cur !~ /.opt[\d]+./){
	    push(@newfiles,$cur);
	}
    }
    @files = @newfiles;


    if($loops > 8)
    {
	print "Stopping recursion at $loops.\n";
    }
}

for(my $index=0; $index<=$#files; $index++){
    open (CKMFFILE, "$files[$index]");

    while(<CKMFFILE> !~ /&flags/){};
    $rt = "";
    $numbf = 0;
    $hfe = "";
    $nw = 0;
    $dt = 0;
    $opt = 0;
    $optl = 0;
    $steps=0;
    $eqsteps=0;
    $searchdata = "";

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

	if($ARGV[0] != "" && $_ =~ m/$ARGV[0]/ && $_ !~ /\#/){
	    $name = $_;
	    chomp $name;
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $val = $line[1];
	    $searchdata .= sprintf "%20s: %30s = %30s\n", "", $name, $val;
	}
	if($_ =~ m/&geometry$/){
	    last;
	}
    }

    $steps_str = "";
    if($steps >= 1000*1000*1000){
	$steps /= int(1000*1000*1000);
	$steps_str = "${steps}G";
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
	$eqsteps_str = "${eqsteps}G";
    } elsif($eqsteps >= 1000*1000){
	$eqsteps /= int(1000*1000);
	$eqsteps_str = "${eqsteps}M";
    } elsif($eqsteps >= 1000){
	$eqsteps /= int(1000);
	$eqsteps_str = "${eqsteps}K";
    } else {
	$eqsteps_str = "$eqsteps";
    }

    if($rt eq "variational"){
	$rt = "v";
    } elsif($rt = "diffusion"){
	$rt = "d";
    }

    printf "%-40s %2s nw=%3i steps=%5s/%-5s dt=%-7s nbf=%-4i opt=%i optl=%i HF= %-20s\n",
    $files[$index],$rt,
    $nw,$eqsteps_str,$steps_str,
    $dt,
    $numbf,$opt,$optl,$hfe;
    print "$searchdata";
    close(CKMFFILE);
}
