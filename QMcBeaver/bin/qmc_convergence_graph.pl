#!/usr/bin/perl

# This script will create a gnuplot graph
# showing the convergence + error bars for a
# set of calculations.
# If the associated .ckmf files are in the same
# directory, it will include some extra info on the graph.
# I programmed this with the latest GNUPLOT on OSX, which
# now has gif output.

use strict;

#hartrees (=1) or kcal/mol (=627.50960803)?
my $units = 627.50960803;

#absolute energies (=0) or relative (=1) to each other?
my $shift = 1;

#keep only 1 line every $drop lines
my $drop = 1;

#should the x axis be samples (=1) or iterations (=0)?
my $xaxis_samples = 0;

my $gnuplot = "/usr/local/bin/gnuplot";
`setenv GDFONTPATH /Library/Fonts:/System/Library/Fonts`;
my $date = `date`;
chomp $date;

if($#ARGV < 0) {
  #if the first arguement is a plotfile.dat, then we'll add new data to it
  #otherwise, we'll create a new plotfile.dat
  print "Usage: qmc_convergence_graph.pl {plotfile.dat || output1.out} [output2.out ...]\n";
  print "Will produce a graph of energy vs num samples or iterations.";
  die;
}

my $y_min = 0;
my $y_max = 0;

my $all_dt = 0;
my $all_nw = 0;

if($ARGV[0] =~ /.dat$/)
{
    print "Adding data to $ARGV[0]\n";
    open (DATFILE, ">>plotfile.dat");
} else {
    print "Truncating plotfile.dat\n";
    open (DATFILE, ">plotfile.dat");
}

my $lastlines = "";
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
    my $line;
    my $more = 1;
    my $eavg;
    my $estd;
    my $iteration;
    my $num_samples;
    my $fordatfile = "";
    my $counter = 0;
    my $numwarnings = 0;
    my $numerrors = 0;
    while(<OUTFILE>){
	if(/WARNING/)
	{
	    if($numwarnings == 0)
	    {
		$base = "?" . $base;
	    }
	    $numwarnings++;
	}
	if(/ERROR/)
	{
	    if($numerrors == 0)
	    {
		$base = "!" . $base;
	    }
	    $numerrors++;
	}
	next if($_ =~ /[a-zA-Z]/ && $_ !~ /Results/);
	chomp;
	my @data = split/[ ]+/;
        #this is the number of data elements per line
	if($#data >= 8){
	    $counter++;
	    $iteration   = $data[1];
	    $eavg        = $data[2];
	    $estd        = $data[3];
	    $num_samples = $data[$#data];
	    next if($num_samples <= 0);

	    #make sure we have the first and last data points included
	    next if($counter%$drop != 0 && $counter != 1 && $iteration%100 == 0);
	    $fordatfile .= sprintf "%20i %20.10f %20.10f %20i\n", $num_samples, $eavg, $estd, $iteration;
	    $line = sprintf "%30s $_\n","$base";
	} elsif(/Results/) {
	    $more = 0;
	}
    }    
    close OUTFILE;
    $lastlines .= "$line";
    my $in_kcal = $eavg*$units;
    printf "%50s %15s %15s E_h=%20.14f E_kcal=%20.10f Err=%i Warn=%i\n","$base","dt=$dt","nw=$nw",$eavg,$in_kcal,$numerrors,$numwarnings;
    #if we are in enhanced text mode, we need to double escape the "_"
    #$base =~ s/_/\\\\_/g;
    printf DATFILE "#%19s %20s %20s %20s\n", "dt=$dt","$base","E=$eavg","";
    print DATFILE "$fordatfile\n\n";
} 
close DATFILE;
print "$lastlines";

#now it's time to generate gnuplot gifs
my @titles;
my $file_name = "qmc_energies.gif";
  
#let's not assume we know what's in the data files
open (DAT_FILE, "plotfile.dat");
my $line = <DAT_FILE>;
while($line){
    chomp $line;
    my @data = split/[= ]+/, $line;
    
    if($line =~ "dt="){
	if($all_dt == -1 && $data[2] != 0){
	    push(@titles,"$data[3], dt=$data[2]"); 
	} else {
	    push(@titles,"$data[3]"); 
	}
    }

    $line = <DAT_FILE>;

    #Make sure we have the last line in a series
    if($line !~ /[0-9]/ && "$data[2]" =~ /[0-9]/){
	if($data[2] < $y_min || $y_min == 0){
	    $y_min = $data[2];
	}
	if($data[2] > $y_max || $y_max == 0){
	    $y_max = $data[2];
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

my $xindex = 4;
my $xlabel = "Num Iterations";
if($xaxis_samples == 1){
    $xindex = 1;
    $xlabel = "Num Samples";
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
set xlabel "$xlabel"
set ylabel "$ylabel"
set grid ytics
set mytics
set ticscale 1.5 0.75
gnuplot_Commands_Done
  
print GNUPLOT "plot ";
for(my $i=0; $i<$#titles; $i++){
    print GNUPLOT " \"plotfile.dat\" index $i using $xindex:(\$2-$shift)*$units:(\$3*$units) title \"$titles[$i]\" with yerrorlines,\\";
    print GNUPLOT "\n";
}
print GNUPLOT " \"plotfile.dat\" index $#titles using $xindex:(\$2-$shift)*$units:(\$3*$units) title \"$titles[$#titles]\" with yerrorlines\n"; 

close (GNUPLOT);
#`/bin/rm $_.dat`;

`open $file_name`;
