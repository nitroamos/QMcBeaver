#!/usr/bin/perl

#assume utilities.pl is in the same directory as summary.pl
my $path = `dirname $0`;
chomp($path);
require "$path/utilities.pl";

my $useVar = 0;

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
#also look at $every in plotter.pl
my $drop = 1;
if($drop != 1){
    print "Keeping only 1 line in $drop\n";
}

my $d = qx! date +%F.%H-%M-%S !;
chomp($d);
#my $gnuplot = "/usr/local/bin/gnuplot";
my $gnuplot = "gnuplot";

#you might need to add this command to your .cshrc
# `setenv GDFONTPATH /Library/Fonts:/System/Library/Fonts`;
# `setenv GDFONTPATH /usr/share/fonts/bitstream-vera/`;

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

getFileList(".out",\@files);

my $y_min = 0;
my $y_max = 0;

my $all_dt = 0;
my $all_nw = 0;

if($ARGV[0] =~ /.dat$/)
{
    print "Adding data to $ARGV[0]\n";
    open (DATFILE, ">>plotfile.dat");
} else {
    #print "Truncating plotfile.dat\n";
    open (DATFILE, ">plotfile.dat");
}

my $num_results;
my $ave_result;

my %dt_ave_results;
my %label;
my %dt_err_results;
my %dt_num;
my %dt_num_results;
my %dt_nme_results;
my %summary;

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

    my $vare = "";

    my $dt = 0;
    my $nw = 0;
    my $effnw = 0;
    my $opt = -1;
    my $isd = -1;
    my $hfe = 0;
    my $numbf = 0;
    my $numci = 0;
    my $use3 = 0;
    open (CKMFFILE, "$base.ckmf");
    while(<CKMFFILE>){
	if($_ =~ /^\#/ && $_ !~ /[A-DF-Za-df-z]+/ && $vare eq ""){
	    chomp;
	    my @line = split/[ ]+/;
	    $vare    = $line[2];
	}

	if($_ =~ m/^\s*run_type\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $isd = $line[1];
	    if($useVar == 0){
		last if($isd eq "variational");
	    }
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
	if($_ =~ m/^\s*use_three_body_jastrow\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $use3 = $line[1];
	}
	if($_ =~ m/&geometry$/){
	    last;
	}
    }

    next if($opt == 1);
    if($useVar == 0){
	next if($isd eq "variational");
    }
    next if($nw != 100);
    next if($dt != 0.01);
    
    while(<CKMFFILE> !~ /Jastrow/){}
    my $numjw = "";
    my $numjwID = 0;
    while(<CKMFFILE>){
	if(/NumberOfParametersOfEachType/){
	    if(!($numjw eq "")){
		$numjw .= ",";
	    }
	    my @line = split/\s+/;
	    my $numthis = "$line[1]";
	    $numjwID += $line[1];
	    for($i=2; $i<=$#line; $i++){
		$numthis .= "$line[$i]";
		$numjwID += $line[$i];
	    }
	    $numjw .= "$numthis";
	}
    }
    $numjw = "$numjwID=$numjw";

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
    my $wallclock = "";
    my $totalclock = "";
    my $effdt = 0;

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
	$wallclock = (split/[ ]+/)[11] if(/Wallclock/);
	$totalclock = (split/[ ]+/)[11] if(/Total Time/);

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
	    $effnw       = $data[4];
	    $effdt       = $data[7];
	    #$num_samples = $data[$#data];
	    $num_samples = $data[8];
	    #printf "$#data $more %20i %20.10f %20.10f %20i\n", $num_samples, $eavg, $estd, $iteration; 
	    next if($num_samples <= 0);

	    #make sure we have the first and last data points included
	    next if($counter%$drop != 0 && $counter != 1 && $iteration%100 == 0);
	    next if($iteration < 0);
	    $fordatfile .= sprintf "%20i %20.10f %20.10f %20i\n", $num_samples, $eavg, $estd, $iteration;
	    $line = sprintf "%30s $_\n","$base";
	} elsif(/Results/) {
	    $more = 0;
	}
    }    
    close OUTFILE;
    $lastlines .= "$line";

    my $key;
    if($vare eq ""){
	$key =  "$hfe&$dt&$numbf&$numjw&$nw&$numci";
    } else {
	$key = "$vare&$dt&$numbf&$numjw&$nw&$numci";
    }

    if($eavg < 0){
	my $weight = $num_samples/100000;

	$dt_ave_results{$key} += $eavg * $weight;
	$dt_num_results{$key} += $weight;
	$dt_num{$key} += 1.0;

	if($estd > 0){
	    $dt_err_results{$key} += $estd * $estd * $weight;
	    $dt_nme_results{$key} += $weight;
	}

	$ave_result += $eavg;
	$num_results++;
    }
    my $in_kcal = $eavg*$units;
    #printf "%50s %15s %15s E_h=%20.14f E_kcal=%20.10f Err=%i Warn=%i\n","$base","dt=$dt","nw=$nw",$eavg,$in_kcal,$numerrors,$numwarnings;

    $summary{$key} .= 
	sprintf "..  %-30s %7s %5s  %-15s %20s %5i %5i %20.10f %4i %5i %10s %10s %10s %15.5f\n",
	"$base","$dt","$effnw","$hfe","$vare",
	$numbf,$numjw,$eavg,
	$numerrors,$numwarnings,$wallclock,$totalclock,
	$iteration,($effdt*$iteration);

    #if we are in enhanced text mode, we need to double escape the "_"
    #$base =~ s/_/\\\\_/g;
    printf DATFILE "#%19s %20s %20s %40s\n", "dt=$dt","$base","E=$eavg","$key";
    print DATFILE "$fordatfile\n\n";
} 
close DATFILE;
print "$lastlines";

foreach $key (sort byenergy keys %dt_ave_results)
{
    if(!exists $label{$key})
    {
	$label{$key} = sprintf "%2i", (scalar keys %label) + 1;
    }
}
printf "ID  %-30s %7s %5s  %-15s %20s %5s %5s %20s %4s %5s %10s %10s %10s %15s\n",
    "File Name","dt","nw","HF Energy",
    "Var. Energy","NumBF",
    "NumJW","Average","Err","Warn",
    "Wall","Total","Iter","effdt*iters";

foreach $sum (sort bydt keys %summary)
{
    $summary{$sum} =~ s/../$label{$sum}/;
    print "$summary{$sum}";
}

if($num_results > 0){
    $ave_result /= $num_results;
    #print "Average result = $ave_result\n";
    
    printf "%5s %10s %20s %5s %5s   %-25s %5s %20s %20s %20s %10s\n",
    "ID",
    "dt","Ref. Energy","Num","NumBF","NumJW","NumW","Average","Error (kcal)","Corr. E.","Weight";
    my %qref;
    my %href;
    my $dtref = 0;
    my $cure = "";
    foreach $key (sort byenergy keys %dt_ave_results)
    {
	my @keydata = split/&/,$key;
	$dt_ave_results{$key} /= $dt_num_results{$key};
	$dt_err_results{$key} = sqrt($dt_err_results{$key}/$dt_nme_results{$key});

	printf "%5i %10s %20s %5i %5i   %-25s %5i %20.10f %20.10f %20.10f %10.5f\n",
	$label{$key},
	"$keydata[1]", "$keydata[0]",
	$dt_num{$key},
	$keydata[2],
	$keydata[3],
	$keydata[4],
	$dt_ave_results{$key},
	$dt_err_results{$key}*$units,
	($keydata[0]-$dt_ave_results{$key}),
	$dt_num_results{$key};
    }

    print "\n\n";
    #matrix output
    #the data is sorted according to dt first
    #we calculate the difference for for all results available
    #but we don't compare calculations if dt and energy are different
    print "DMC comparisons:\n";
    foreach $row (sort bydt keys %dt_ave_results)
    {
	my $newrow = 1;
	my @rowdata = split/&/,$row;
	my $rowJW = (split/=/,$rowdata[3])[0];
	foreach $col (sort bydt keys %dt_ave_results)
	{
	    my @coldata = split/&/,$col;
	    my $colJW = (split/=/,$coldata[3])[0];

	    ($rMult,$cMult) = getFormula($row,$col);

	    my $r = $dt_ave_results{$row};
	    my $c = $dt_ave_results{$col};

	    #the results are not comparable if either is zero
	    next if($rMult == 0 || $cMult == 0);
	    #otherwise we'll get two of every comparison
	    next if($r < $c);

	    my $diff = $r*$rMult - $c*$cMult;
	    my $diffe = $dt_err_results{$row}*$rMult +
		$dt_err_results{$col}*$cMult;
	    
	    $diff *= -$units;
	    $diffe *= $units;
	    
	    if($newrow == 1){
		printf "%3i) %15.10f x $rMult %9s %5s %5s %2i | ",$label{$row},$r,$rowdata[1],$rowdata[2],$rowJW,$rowdata[5];
		$newrow = 0;
	    } else {
		printf "%49s | "," ";
	    }
	    printf "%3i) %15.10f x $cMult %9s %5s %5s %2i | ",$label{$col},$c,$coldata[1],$coldata[2],$colJW,$coldata[5];
	    printf " %9.3f",$diff;
	    printf " +/- %-9.3f", $diffe;
	    printf "\n";	    
	}
    }
    print "\n\n";
}
