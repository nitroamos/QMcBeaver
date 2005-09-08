#!/usr/bin/perl
use strict;

# This script, a work in progress, I eventually intend to use in some kind
# of regression testing set up. Essentially, this creates an easy way to
# run several different input files (specifed as input to this script)
# on all the executables in $bindir.
#
# The usefulness is that it allows you to easily monitor changes in
# energy and time taken over the different executables.
#
# With Chip's changes to the exe naming style, I have not figured out how
# to get this script to identify which exes that it can run.

my $bindir = "../bin/";
my $os = "";
my $d = qx! date !;
my $p = qx! pwd !;
my @exelist;

if($os == ""){
  if($ENV{OSTYPE} =~ /irix/){
    #$os = "sgi";
    $os = "x";
  } elsif($ENV{OSTYPE} =~ /linux/){
    #$os = "linux";
    $os = "x";
  } else {
    #print "unknown computer, please specify exe suffix\n";
    #die;
    $os = "x";
  }
}

if($#ARGV < 0) {
    print "Usage: runallqmcbeavers.pl input1.ckmf [input2.ckmf ...]\n";
    die;
}

my @exelist = qx! ls $bindir !;
chomp @exelist;

system "/bin/rm runall_results.txt";
open (RESULTS, '>runall_results.txt');
print RESULTS "Results compiled on $d$p\n";

foreach my $inputfile (@ARGV){
  my @inputtitle = split/\./,$inputfile;
  if($inputtitle[$#inputtitle] != "ckmf") {next;}

  foreach my $exe (@exelist){
    if(-d $exe){ next; }
    my @exetitle = split/\./,$exe;
    if($exetitle[$#exetitle] !~ /$os/) { next; }

    print "Running $bindir$exe $inputfile...\n";
    `$bindir$exe $inputfile`;

    my $rsltsfile = $inputfile;
    $rsltsfile =~ s/ckmf/rslts/;
    if(!(-e $rsltsfile)){
      print "$bindir$exe $inputfile did not complete\n";
    }

    open (OUTPUT, "$rsltsfile");

    my $energy;
    my $time;
    my $line = <OUTPUT>;
    while(defined $line && ($line !~ "Energy" || $line !~ "Total Time")){
      if($line =~ "- Energy -"){
        $line = <OUTPUT>;
        $energy = $line;
      } elsif($line =~ "Total Time"){
        my @temp = split/ +/,$line;
        $time = $temp[2];
        $line = <OUTPUT>;
      } else {
        $line = <OUTPUT>;
      }
    }

    print RESULTS "$bindir$exe $inputfile gave the following results:\n";
    print RESULTS "$inputfile $exetitle[$#exetitle - 1] Time Taken: $time ms\n";
    print RESULTS "$inputfile $exetitle[$#exetitle - 1] Final Energy: $energy\n";

    close OUTPUT;
  }
}

close RESULTS;
