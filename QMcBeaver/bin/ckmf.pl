#!/usr/bin/perl
#this script will help you examine many ckmf files simultaneously

my $path = `dirname $0`;
chomp($path);
require "$path/utilities.pl";

my $includeRestarts = 0;

my @files = @ARGV;
if($#files < 0)
{
    push(@files,".");
}
getFileList(".ckmf",\@files);

my $curTime = qx! date +%s !;

printf "%-30s %2s %3s %11s %-7s %-4s %4s %4s  %-15s %8s %8s\n",
    "Name","  ","NW","EQ/Steps","dt","nbf","opt","optl","HF","Age","Size";
my $base = "";
for(my $index=0; $index<=$#files; $index++){
    if(!$includeRestarts){
	next if($files[$index] =~ /.\d\d.ckmf$/);
    }
    $base = substr($files[$index],0,-5);
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
    $iseed=0;
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
	if($_ =~ m/^\s*iseed\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $iseed = $line[1];
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

    my $outModTime = "";
    my $outSize    = "";
    my $ovData     = "";
    my $failed     = 0;

    if(-e "$base.out"){
	my $data = `/bin/ls -lh --time-style=+%s $base.out`;
	my @list = split/ +/,$data;

	$data = `grep failed $base.out`;
	$failed = 1 if(length($data) > 0);

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
    }

    printf "%-30s %2s %3i %5s/%-5s %-7s %-4i %4i %4i %-15s",
    $base,$rt,
    $nw,$eqsteps_str,$steps_str,
    $dt,
    $numbf,$opt,$optl,$hfe;
    printf " %10s", $outModTime;
    printf " %7s", $outSize;
    printf " %40s", $ovData;
    if($iseed != 0){
	print " iseed = $iseed";
    }
    print "\n";

    print "$searchdata";
    close(CKMFFILE);
}
