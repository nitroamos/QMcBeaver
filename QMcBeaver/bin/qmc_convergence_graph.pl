#!/usr/bin/perl

# This script will create a gnuplot graph
# showing the convergence + error bars for a
# set of calculations.
# If the associated .ckmf files are in the same
# directory, it will include some extra info on the graph.
# I programmed this with the latest GNUPLOT on OSX, which
# now has gif output.

#use strict;

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

my @files = @ARGV;
if($#ARGV < 0)
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
	if(-d $cur && $cur !~ /src$/ && $cur !~ /bin$/ && $cur !~ /include$/)
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
	} elsif($cur =~ /.out$/){
	    push(@newfiles,$cur);
	}
    }
    @files = @newfiles;


    if($loops > 8)
    {
	print "Stopping recursion at $loops.\n";
    }
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

my $num_results;
my $ave_result;

my %dt_ave_results;
my %dt_err_results;
my %dt_num_results;

printf "%50s %7s %7s %20s %20s %10s %10s\n","File Name","dt","nw","HF Energy","Average","Errors","Warnings";
my $lastlines = "";
for(my $index=0; $index<=$#files; $index++){
    my $cur = $files[$index];
    next if(!(-f $cur));
    my $base = "";
    if($cur =~ /.out$/){
	$base = substr($cur,0,-4);
    } elsif($cur =~ /.qmc$/){
	$base = substr($cur,0,-4);
    } else {
	next;
    }
    my $dt = 0;
    my $nw = 0;
    my $opt = -1;
    my $isd = -1;
    my $hfe = 0;
    open (CKMFFILE, "$base.ckmf");
    while(<CKMFFILE>){
	if($_ =~ m/^\s*run_type\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $isd = $line[1];
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
	if($_ =~ m/^\s*energy\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $hfe = $line[1];
	}
	if($_ =~ m/&geometry$/){
	    last;
	}
    }

    #next if($dt ne "0.001");
    next if($opt == 1);
    next if($isd eq "variational");

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

    open (OUTFILE,  "$cur");
    my $line;
    my @data;
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
		#$base = "?" . $base;
	    }
	    $numwarnings++;
	}
	if(/ERROR/)
	{
	    if($numerrors == 0)
	    {
		#$base = "!" . $base;
	    }
	    $numerrors++;
	}

	next if($_ =~ /[=]/ && $_ !~ /Results/);
	chomp;
	@data = split/[ ]+/;
        #this is the number of data elements per line
	if($#data >= 8 && $more){
	    $counter++;
	    $iteration   = $data[1];
	    $eavg        = $data[2];
	    $estd        = $data[3];
	    $num_samples = $data[$#data];
	    #printf "$#data $more %20i %20.10f %20.10f %20i\n", $num_samples, $eavg, $estd, $iteration; 
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
    if($eavg < 0){
	my $key = "$hfe&$dt";
	$dt_ave_results{$key} += $eavg;
	$dt_err_results{$key} += $estd;
	$dt_num_results{$key} += 1;
	$ave_result += $eavg;
	$num_results++;
    }
    my $in_kcal = $eavg*$units;
    #printf "%50s %15s %15s E_h=%20.14f E_kcal=%20.10f Err=%i Warn=%i\n","$base","dt=$dt","nw=$nw",$eavg,$in_kcal,$numerrors,$numwarnings;
    printf "%50s %7s %7s %20s %20.10f %10i %10i\n","$base","$dt","$nw","$hfe",$eavg,$numerrors,$numwarnings;
    #if we are in enhanced text mode, we need to double escape the "_"
    #$base =~ s/_/\\\\_/g;
    printf DATFILE "#%19s %20s %20s %20s\n", "dt=$dt","$base","E=$eavg","";
    print DATFILE "$fordatfile\n\n";
} 
close DATFILE;
print "$lastlines";

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
	$bdata[1] <=> $adata[1];
    } else {
	$bdata[0] <=> $adata[0];
    }
}

if($num_results > 0){
    $ave_result /= $num_results;
    #print "Average result = $ave_result\n";
    
    printf "%10s %20s %20s %20s %20s %20s %20s %20s %10s\n",
    "dt","HF Energy","Average","Error (kcal)","dt Diff. (kcal)","Diff. (kcal)","HF Diff. (kcal)","Corr. E.","Number";
    my %qref;
    my %href;
    my $dtref = 0;
    my $cure = "";
    foreach $key (sort byenergy keys %dt_ave_results)
    {
	my @keydata = split/&/,$key;
	$dt_ave_results{$key} /= $dt_num_results{$key};
	$dt_err_results{$key} /= $dt_num_results{$key};

	if(! exists $qref{$keydata[1]}){
	    $qref{$keydata[1]} = $dt_ave_results{$key};
	    $href{$keydata[1]} = $keydata[0];
	    #print "Setting $keydata[1] reference to $qref{$keydata[1]}\n";
	} elsif(abs($qref{$keydata[1]} - $dt_ave_results{$key}) > 5){
	    #if the energies are different by more than 10 hartrees, assume we want to reset
	    $ratio = $dt_ave_results{$key} / $qref{$keydata[1]};
	    $intr = int($ratio);
	    $intdiff = $ratio - $intr;
	    #print "int = $intr diff = $intdiff\n";
	    if($intr != 1 && $intdiff < 0.01)
	    {
		$qref{$keydata[1]} *= $intr;
		$href{$keydata[1]} *= $intr;
	    } else {
		#assume it's completely different
		$qref{$keydata[1]} = $dt_ave_results{$key};
		$href{$keydata[1]} = $keydata[0];
		#print "Setting $keydata[1] reference to $qref{$keydata[1]}\n";
	    }
	}

	if($cure ne $keydata[0]){
	    $cure = $keydata[0];
	    $dtref = $dt_ave_results{$key};
	}
	#print "ref = $qref{$keydata[1]}\n";

	printf "%10s %20s %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %10i\n",
	"$keydata[1]", "$keydata[0]",
	$dt_ave_results{$key},
	$dt_err_results{$key}*$units,
	($dtref - $dt_ave_results{$key})*$units,
	($qref{$keydata[1]} - $dt_ave_results{$key})*$units,
	($href{$keydata[1]} - $keydata[0])*$units,
	($keydata[0]-$dt_ave_results{$key}),
	$dt_num_results{$key};
    }

    #matrix output
    #the data is sorted according to dt first
    #we calculate the difference for for all results available
    #but we don't compare calculations if dt and energy are different
    print "\n";
    foreach $col (sort bydt keys %dt_ave_results)
    {
	my @coldata = split/&/,$col;
	printf "%9.4f   ",$coldata[0];
    }
    print "\n";
    foreach $col (sort bydt keys %dt_ave_results)
    {
	my @coldata = split/&/,$col;
	printf "%9s   ",$coldata[1];
    }
    print "\n";
    foreach $col (sort bydt keys %dt_ave_results)
    {
	printf "%9s-- ","---------";
    }
    print "\n";
    foreach $row (sort bydt keys %dt_ave_results)
    {
	my @rowdata = split/&/,$row;
	foreach $col (sort bydt keys %dt_ave_results)
	{
	    my @coldata = split/&/,$col;
	    my $print = 1;
	    if($rowdata[0] != $coldata[0]){
		$print = 0 if($rowdata[0] <= $coldata[0] ||
			      $rowdata[1] != $coldata[1]);
	    } else {
		$print = 0 if($rowdata[1] <= $coldata[1]);
	    }

	    my $r = $dt_ave_results{$row};
	    my $c = $dt_ave_results{$col};
	    my $diff = $r - $c;

	    my $intr = 0;
	    if(abs($r/$c) > 1.5)
	    {
		#if the energies are different by more than 10 hartrees, assume we want to reset
		$ratio = $r / $c;
		$intr = int($ratio);
		$intdiff = $ratio - $intr;

		if($intr != 1 && $intdiff < 0.1)
		{
		    $diff = $r - $c * $intr;
		}
	    } elsif($c/$r > 1.5) {
		#if the energies are different by more than 10 hartrees, assume we want to reset
		$ratio = $c / $r;
		$intr = int($ratio);
		$intdiff = $ratio - $intr;

		if($intr != 1 && $intdiff < 0.1)
		{
		    $diff = $intr * $r - $c;
		}
	    }
	    
	    $diff *= $units;

	    if($print)
	    {
		if(abs($diff) > 1e3) {
		    printf "%8.2fK   ",($diff/1000.0);
		} else {
		    if($intr > 1)
		    {
			printf "%9.3f %i ",$diff,$intr;
		    } else {
			printf "%9.3f   ",$diff;
		    }
		}
	    } else {
		printf "%9s   "," ";
	    }
	}
	printf "| %9.4f %9s\n",$rowdata[0],$rowdata[1];
    }
    print "\n\n";
}

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
