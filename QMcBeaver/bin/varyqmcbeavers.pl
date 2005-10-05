#!/usr/bin/perl
use strict;

# This script, a work in progress, is meant to systematically vary any QMcBeaver ckmf parameter
# run, and collect results out of the stream dumped to stdout.
#
# The stdout data of interest needs to be marked with a '#' (or whatever is in $oe_exe that doesn't
# conflict with regex). This script will scan through the output grabbing $numout data points marked with
# $find_tag.
#
# Put the name of the parameter to be varied in $param_name.
# You'll need to modify this code to indicate how you want that parameter to be varied.

#for a smoothly varying parameter
my $d1 = 1;
my $delta = 1;
my $numruns = 2;
my $numout = 2;
#my $param_name = "walkers_per_pass";
#my $param_name2 = "number_of_walkers";
my $param_name = "iseed";
my $do_random = 1;
my $delete_files = 1;

my $find_tag = "#";

my $bindir = "~/QMcBeaver/bin/";
my $os_exe = "";
my $base_input = "";
my $base_base = "";
my $results_file = "";
my $d = qx! date !;
my $time_start = qx! date +%s !;
my $p = qx! pwd !;
my @exelist;
my %inputlist;

if($os_exe == ""){
  if($ENV{OSTYPE} =~ /irix/){
    #$os_exe = "sgi";
    $os_exe = "x";
  } elsif($ENV{OSTYPE} =~ /linux/){
    #$os_exe = "linux";
    $os_exe = "x";
  } else {
    #print "unknown computer, please specify exe suffix\n";
    #die;
    $os_exe = "exe";
  }
}

my @exelist;
if($#ARGV < 0) {
	print "usage: <script> <input>.ckmf [exe list].x\n";
} elsif($#ARGV == 0){
	$base_input = $ARGV[0];
	@exelist = qx! ls $bindir !;
	chomp @exelist;
} else {
	$base_input = $ARGV[0];
	@exelist = @ARGV[1 .. $#ARGV];
}

my @temp = split/\./,$base_input;
$base_base = join("",@temp[0 .. $#temp-1]);

$results_file = "runall.$base_base.$param_name.txt";

if(-e $results_file){
  my @templist = qx! ls runall.$base_base.$param_name* !;
  my $index = $#templist + 1;
	$results_file = "runall.$base_base.$param_name.$index.txt";
}

open (RESULTS, ">$results_file");
print RESULTS "Results compiled on $d$p";
print RESULTS "Input: $base_input\n";
print RESULTS "Varied: $param_name\n";

for(my $i=0; $i<$numruns; $i++){
  my $value;
  if($do_random){
    $value = int(rand(-2140000000));
	} else {
  	$value = $d1 + $i*$delta;
  }
  my $inputfile = "${base_base}_$value.ckmf";
	
	if(-e $inputfile){
		#next;
		system "/bin/rm -f $inputfile";
	}
	
	open(BASE_CKMF,"$base_input");
	open(INPUT,">$inputfile");
		
	my $line = <BASE_CKMF>;
	while(defined $line){
		if($line =~ /$param_name/){
			print INPUT $line;
			$line = <BASE_CKMF>;
			print INPUT " $value\n";
			$line = <BASE_CKMF>;
		} else {
			print INPUT $line;
			$line = <BASE_CKMF>;
		}
	}

	close BASE_CKMF;
	close INPUT;
	
	$inputlist{ $value } = $inputfile;
}

foreach my $exe (@exelist){
	if(-d $exe){ next; }
	my @exetitle = split/\./,$exe;
	if($exetitle[$#exetitle] !~ /^$os_exe/) { next; }

	print "Running $bindir$exe\n";
	print RESULTS "\nResults from $bindir$exe\n";

	foreach my $value (sort { $a <=> $b } keys %inputlist) {
		my $outCounter = 0;
		my $input = $inputlist{$value};
		my $rsltsfile = $input;
		$rsltsfile =~ s/ckmf/out/;
		
    print "Doing $value\n";
		system "/bin/rm -f $rsltsfile";
		system "$bindir$exe $input > $rsltsfile";
		
		if(!(-e $rsltsfile)){
			print "$bindir$exe $input did not complete\n";
		}

		open (OUTPUT, "$rsltsfile");

		print RESULTS "$value ";
		my $line = <OUTPUT>;
		while(defined $line){
			my @pieces = split/ +/,$line;
			chomp(@pieces);
			for(my $i=0; $i<$#pieces; $i++){
				if($pieces[$i] =~ /$find_tag/ && $outCounter < $numout){
					print RESULTS "$pieces[$i+1] ";
					$outCounter++;
				}			
			}
			$line = <OUTPUT>;
		}

		print RESULTS "\n";	
		close OUTPUT;
    `rm -f ${base_base}_$value.out`;
    `rm -f ${base_base}_$value.qmc`;
    `rm -f ${base_base}_$value.rslts`;
    `rm -f ${base_base}_$value.energy.0`;
    `rm -f ${base_base}_$value.01.ckmf`;
  }
  `rm -f shader_*`;
}
my $time_stop = qx! date +%s !;
my $total_time = $time_stop - $time_start;
print RESULTS "Time Taken: $total_time seconds\n";
print "Time Taken: $total_time seconds\n";
close RESULTS;

if($delete_files){
  foreach my $value (keys %inputlist) {
    `rm -f ${base_base}_$value.*`;
  }
}