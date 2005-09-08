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

my $bindir = "~/QMcBeaver/bin/";
my $os_exe = "";
my $base_input = "";
my $base_base = "";
my $d = qx! date !;
my $p = qx! pwd !;
my @exelist;
my %inputlist;

#for a smoothly varying parameter
my $d1 = 5;
my $delta = 5;
my $numruns = 40;
my $numout = 2;
#my $param_name = "programmers_longs";
my $param_name = "iseed";
my $find_tag = "#";

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
	$os_exe = "x";
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

system "/bin/rm -f runall_results.txt";
open (RESULTS, '>runall_results.txt');
print RESULTS "Results compiled on $d$p\n";

for(my $i=0; $i<$numruns; $i++){
	#my $value = $d1 + $i*$delta;
	my $value = int(rand(-2140000000));
	my $inputfile = "${base_base}.$value.ckmf";
	
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

	print "Running $bindir$exe...\n";
	print RESULTS "\nResults from $bindir$exe...\n";

	
	foreach my $value (sort { $a <=> $b } keys %inputlist) {
		my $outCounter = 0;
		my $input = $inputlist{$value};
		my $rsltsfile = $input;
		$rsltsfile =~ s/ckmf/out/;
		
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
    }
}

close RESULTS;
