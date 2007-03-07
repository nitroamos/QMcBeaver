#!/usr/bin/perl

# This script will create a gnuplot graph
# showing the convergence + error bars for a
# set of calculations.
# If the associated .ckmf files are in the same
# directory, it will include some extra info on the graph.
# I programmed this with the latest GNUPLOT on OSX, which
# now has gif output.

use strict;

#hartrees or kcal/mol?
my $units = 627.50960803;

#absolute energies or relative to each other?
my $shift = 1;

#keep only 1 line every $drop lines
my $drop = 1;

my $gnuplot = "/usr/local/bin/gnuplot";
`setenv GDFONTPATH /Library/Fonts:/System/Library/Fonts`;
my $date = `date`;
chomp $date;

if($#ARGV < 0) {
  print "Usage: qmc_convergence_graph.pl output1.out [output2.out ...]\n";
  print "Will produce a graph of energy vs num samples";
  die;
}

my $y_min = 0;
my $y_max = 0;

my $all_dt = 0;
my $all_nw = 0;

open (DATFILE, ">plotfile.dat");
for(my $index=0; $index<=$#ARGV; $index++){
    if(!(-f $ARGV[$index]) || $ARGV[$index] !~ /.out$/){ next; }
    my $base = substr($ARGV[$index],0,-4);
    
    my $dt = 0;
    my $nw = 0;
    open (CKMFFILE, "$base.ckmf");
    while(<CKMFFILE>){
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
	if($nw != 0 && $dt != 0){
	    last;
	}
    }

    if($all_dt == 0){
	$all_dt = $dt;
    } elsif($all_dt != $dt){
	$all_dt = -1;
    }
    if($all_nw == 0){
	$all_nw = $nw;
    } elsif($all_nw != $nw){
	$all_nw = -1;
    }
    close CKMFFILE;

    open (OUTFILE,  "$base.out");
    my $line = <OUTFILE>;
    while($line !~ "Iteration"){
	$line = <OUTFILE>;
    }

    my $more = 1;
    my $eavg;
    my $estd;
    my $iteration;
    my $num_samples;
    my $fordatfile = "";
    my $counter = 0;
    my $warned = 0;
    while($line and $more == 1){
	$line = <OUTFILE>;
	if($line =~ /WARNING/ && $warned == 0){
	    $base = "*" . $base;
	    $warned = 1;
	}
	chomp $line;
	my @data = split/[ ]+/,$line;
	#print "line $#data is $line\n";
	if($#data == 8){
	    $counter++;
	    $iteration   = $data[1];
	    $eavg        = $data[2];
	    $estd        = $data[3];
	    $num_samples = $data[8];
	    #make sure we have the first and last data points included
	    next if($counter%$drop != 0 && $counter != 1 && $iteration%100 == 0);
	    $fordatfile .= sprintf "%20i %20.10f %20.10f %20i\n", $num_samples, $eavg, $estd, $iteration;
	} elsif($line =~ /Results/) {
	    $more = 0;
	}
    }    
    close OUTFILE;

    my $in_kcal = $eavg*$units;
    printf "%50s %15s %15s E_h=%20.14f E_kcal=%20.10f\n","$base","dt=$dt","nw=$nw",$eavg,$in_kcal;
    #if we are in enhanced text mode, we need to double escape the "_"
    #$base =~ s/_/\\\\_/g;
    printf DATFILE "#%19s %20s %20s %20s\n", "dt=$dt","$base","E=$eavg","";
    print DATFILE "$fordatfile\n\n";

    if($eavg < $y_min || $y_min == 0){
	$y_min = $eavg;
    }
    if($eavg > $y_max || $y_max == 0){
	$y_max = $eavg;
    }
} 
close DATFILE;

#now it's time to generate gnuplot gifs
my @titles;
my $file_name = "qmc_energies.gif";
  
#let's not assume we know what's in the data files
open (DAT_FILE, "plotfile.dat");
while(my $line = <DAT_FILE>){
    if($line =~ "dt="){
	chomp $line;
	my @data = split/[= ]+/, $line;
	if($all_dt == -1 && $data[2] != 0){
	    push(@titles,"$data[3], dt=$data[2]"); 
	} else {
	    push(@titles,"$data[3]"); 
	}
    }
}

if($shift == 1){
    $shift = $y_min;
    $y_max = ($y_max-$y_min)*$units;
    $y_min = 0;
} else {
    $shift = 0
}

my $ylabel = "Energy";
if($units == 1){
    $ylabel .= " (au)";
} else {
    $ylabel .= " (kcal/mol)";
}
my $space = 0.1*($y_min - $y_max);
$y_min += $space;
$y_max -= $space;

my $title_extra = "";
if($all_dt != -1){
    $title_extra .= ", dt=$all_dt";
}
if($all_nw != -1){
    $title_extra .= ", nw=$all_nw";
}

`/bin/rm -f $file_name`;
open(GNUPLOT, "|$gnuplot");
print GNUPLOT <<gnuplot_Commands_Done;
#fonts with extensions "ttf" and "dfont" will work
#here is a list of available fonts: Chalkboard Helvetica Times
#Courier Monaco LucidaGrande
set term gif crop enhanced font 'Monaco' 8
set output "$file_name"

#set term png
#set terminal png medium
set size 0.9,1

set nokey
set key outside below box noenhanced
set yrange[$y_min:$y_max]
set title "QMC Runs${title_extra}\\n{/=8${date}}"
set xlabel "Num Samples"
set ylabel "$ylabel"
set grid ytics
set mytics
set ticscale 1.5 0.75
gnuplot_Commands_Done
  
print GNUPLOT "plot ";
for(my $i=0; $i<$#titles; $i++){
    print GNUPLOT " \"plotfile.dat\" index $i using 1:(\$2-$shift)*$units:(\$3*$units) title \"$titles[$i]\" with yerrorlines,\\";
    print GNUPLOT "\n";
}
print GNUPLOT " \"plotfile.dat\" index $#titles using 1:(\$2-$shift)*$units:(\$3*$units) title \"$titles[$#titles]\" with yerrorlines\n"; 

close (GNUPLOT);
#`/bin/rm $_.dat`;

`open $file_name`;
