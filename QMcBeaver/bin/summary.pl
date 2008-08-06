#!/usr/bin/perl

#assume utilities.pl is in the same directory as summary.pl
my $path = `dirname $0`;
chomp($path);
require "$path/utilities.pl";

my $useVar     = 0;
my $dtFilter   = 0.01;
my $orbFilter  = 1;
my $compareE   = 0;
my $sumResults = 1;
my $latexHelp  = 0;
my @fileFilters;
my @exclusionFilters;

#my $extraTag  = "trail_eps2";

# This script will create a gnuplot graph
# showing the convergence + error bars for a
# set of calculations.
# If the associated .ckmf files are in the same
# directory, it will include some extra info on the graph.
# I programmed this with the latest GNUPLOT on OSX, which
# now has gif output.

#use strict;

#hartrees (=1) or kcal/mol (=627.50960803) or eV (27.211399)?
#my $units = 27.211399;
my $units = 627.50960803;
my $unitsL = "kcal/mol";


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

while($#ARGV >= 0 && $ARGV[0] =~ /^-/){
    $type = shift(@ARGV);
    $param = "";

    if($type eq "-v"){
	$param = shift(@ARGV);
	$useVar = $param if($param == 1 || $param == 0);
	print "Using useVar = $useVar  $useVar\n";
    } elsif($type eq "-t"){
	$param = shift(@ARGV);
	$dtFilter = $param if($param >= 0);
	print "Using dt filter $dtFilter\n";
    } elsif($type eq "-f"){
	$param = shift(@ARGV);
	push(@fileFilters,$param);
	print "Adding file filter $param\n";
    } elsif($type eq "-x"){
	$param = shift(@ARGV);
	push(@exclusionFilters,$param);
	print "Adding file exclusion filter $param\n";
    } elsif($type eq "-u"){
	$param = shift(@ARGV);
	if($param == 0){
	    $units  = 627.50960803;
	    $unitsL = "kcal/mol";
	} elsif($param == 1) {
	    $units  = 27.211399;
	    $unitsL = "eV"; 
	} elsif($param == 2) {
	    $units  = 2625.5002;
	    $unitsL = "kJ/mol"; 
	} elsif($param == 3) {
	    $units  = 219474.63;
	    $unitsL = "cm^-1"; 
	} elsif($param == 4) {
	    $units  = 1;
	    $unitsL = "au"; 
	}
	print "Using $unitsL energy units, conversion = $units\n";
    } elsif($type eq "-o"){
	$orbFilter = 0;
	if($orbFilter == 1){
	    print "Filtering to only include balanced orbitals\n";
	} else {
	    print "Not filtering results based on orbital usage.\n";
	}
    } elsif($type eq "-c"){
	$compareE = 1;
	print "Comparing with reference energies.\n";
    } elsif($type eq "-s"){
	$sumResults = 0;
	print "Extended report.\n";
    } elsif($type eq "-l"){
	$latexHelp = 1;
	print "LaTex Helper.\n";
    } else {
	print "Unrecognized option: $type\n";
	die;
    }
}

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

#getFileList(".out",\@files);
getFileList(".qmc",\@files);

if($ARGV[0] =~ /.dat$/)
{
    print "Adding data to $ARGV[0]\n";
    open (DATFILE, ">>plotfile.dat");
} else {
    #print "Truncating plotfile.dat\n";
    open (DATFILE, ">plotfile.dat");
}

my $Cnormal = "\x1b[0m";
my $Chilite = "\x1b[37m";

my $lenLong = 0;
my $num_results;
my $ave_result;
my $headerLine = "";

my %dt_ave_results;
my %label;
my %dt_err_results;
my %dt_num;
my %dt_num_results;
my %dt_nme_results;
my %summary;
my %shortnames;
my %referenceE = ();

my $lastlines = "";
for(my $index=0; $index<=$#files; $index++){
    my $cur = $files[$index];
    next if(!(-f $cur));

    my $isIncluded = 1;
    my $filterMatch = 0;
    foreach $filter (@fileFilters){
	#we only include a file if it matches one of the filters
	$filterMatch = 1  if($cur =~ /$filter/);
    }
    $isIncluded = 0 if($filterMatch == 0 && $#fileFilters >= 0);

    foreach $filter (@exclusionFilters){
	#exclude a file if it matches one of the exclusion filters
	#print "filter = $filter cur = $cur\n";
	$isIncluded = 0 if($cur =~ /$filter/);
    }
    next if($isIncluded == 0);

    my $base = "";
    if($cur =~ /.out$/){
	$base = substr($cur,0,-4);
    } elsif($cur =~ /.qmc$/){
	$base = substr($cur,0,-4);
    } else {
	next;
    }
    my $short = `basename $base`;
    chomp($short);
    #remove the restart index
    $short = $1 if($short =~ /([\w\d]+)\.\d\d$/);
    #remove any _\d at the end
    $short = $1 if($short =~ /([\w\d]+)_[\d]+$/);

    my $vare = "";

    my $dt = 0;
    my $nw = 0;
    my $effnw = 0;
    my $opt = -1;
    my $isd = -1;
    my $hfe = 0;
    my $numbf = 0;
    my $numci = 0;
    my $numor = 0;
    my $use3 = 0;
    my $extraVal = 0;
    my %refEnergies = ();

    open (CKMFFILE, "$base.ckmf");
    while(<CKMFFILE>){
	if($_ =~ /^\#/ && $_ !~ /[A-DF-Za-df-z]+/ && $vare eq ""){
	    chomp;
	    my @line = split/[ ]+/;
	    #This is from the header; the top energy is the best
	    $vare    = $line[2];
	    $refEnergies{"VMC"} = $vare;
	}

	$refEnergies{"RHF"} = (split/\s+/)[5] if(/FINAL RHF ENERGY/);
	$refEnergies{"RHF"} = (split/\s+/)[5] if(/FINAL ROHF ENERGY/);
	$refEnergies{"GVB"} = (split/\s+/)[5] if(/FINAL GVB ENERGY/);
	$refEnergies{"CI"} = (split/\s+/)[4] if(/^\#/ && /STATE/ && /ENERGY/);
	$refEnergies{"$1"} = (split/\s+/)[3] if(/^\#/ && /\s+([\w\d\(\)]+)\s+ENERGY:/);

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
	if($_ =~ m/^\s*norbitals\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $numor = $line[1];
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
	if($extraTag ne ""){
	    if($_ =~ m/^\s*$extraTag\s*$/){
		$_ = <CKMFFILE>;
		chomp;
		my @line = split/[ ]+/;
		$extraVal = $line[1];
	    }
	}
	if($_ =~ m/&geometry$/){
	    last;
	}
    }

    next if($opt == 1);
    if($useVar == 0){
	next if($isd eq "variational");
    }
    if($nw < 100){
	print "Not including $base because it has $nw walkers.\n";
	next;
    }
    next if($dtFilter != 0 && $dt != $dtFilter);
    
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
    close CKMFFILE;

    open(RUNFILE, "$base.run");
    my $machine = "";
    while(<RUNFILE>){
	if(/lamboot/){
	    $machine = "m";
	} elsif(/machinefile/){
	    $machine = "h";
	}
    }
    close(RUNFILE);

    open (QMCFILE,  "$cur");
    my $line;
    my @data;
    my $more = 1;
    my $eavg;
    my $estd;
    my $iteration;
    my $num_samples = 0.00001;
    my $fordatfile = "";
    my $counter = 0;
    my $wallclock = "";
    my $totalclock = "";
    my $sampleclock = "";
    my $effdt = 0;
    my $sampleVar = 0;
    my $sampleVarCorLen = 0;
    my $corLength = 0;

    while(<QMCFILE>){
	$headerLine = $_ if(/iteration/ && /Eavg/ && /Samples/);

	#this is to avoid processing lines with warnings
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

	    #In the old format, this was trial energy
	    #In the new format, this is the acceptance probability, which uses parenthesis
	    if($data[6] =~ /\(/)
	    {
		#new output format
		$num_samples = $data[4];
		$corLength   = $data[5];
		$effnw       = $data[7];
		if($isd eq "variational"){
		    $effdt       = $dt;
		} else {
		    $effdt       = $data[10];
		}
	    } else {
		$effnw       = $data[4];
		#old output format
		if($isd eq "variational"){
		    $effdt       = $data[5];
		    $num_samples = $data[6];
		} else {
		    $effdt       = $data[7];
		    $num_samples = $data[8];
		    
		}
	    }
	    
	    #this is equal to sample variance * correlation length
	    $sampleVarCorLen = $estd * $estd * $num_samples;

	    next if($num_samples <= 0);

	    #make sure we have the first and last data points included
	    next if($counter%$drop != 0 && $counter != 1 && $iteration%100 == 0);
	    next if($iteration < 0);
	    $fordatfile .= sprintf "%20i %20.10f %20.10f %20i\n", $num_samples, $eavg, $estd, $iteration;
	    if($extraTag ne ""){
		$line = sprintf "%30s $_ %15s\n","$base","$extraVal";
	    } else {
		$line = sprintf "%30s $_\n","$base";
	    }
	    
	} elsif(/Results/) {
	    $more = 0;
	}
    }    
    close QMCFILE;

    my @times = `grep Time $base.out`;
    $wallclock = (split/\s+/,$times[0])[11];
    $totalclock = (split/\s+/,$times[1])[11];
    $sampleclock = (split/\s+/,`grep "per sample per" $base.out`)[8];
    chomp($sampleclock);
    my $numwarnings = `grep WARNING $base.out | wc -l`;
    #This is "number of warnings per 1000 samples"
    $numwarnings /= $num_samples/1e3;
    $numwarnings = sprintf "%4.2f", $numwarnings;
    $numwarnings = " $Chilite$numwarnings$Cnormal" if($numwarnings > 0.5);
    my $numerrors = `grep ERROR $base.out | wc -l`;
    $numerrors /= $num_samples/1e3;
    $numerrors = sprintf "%4.2f",$numerrors;

    $lastlines .= "$line";
    my $key;
    if($vare eq ""){
	#use the value for energy in the key
	$key =  "$hfe&$dt&$numbf&$numjw&$nw&$numci&$numor";
    } else {
	#use the variational energy from the header in the key
	$key = "$vare&$dt&$numbf&$numjw&$nw&$numci&$numor";
    }

    if(exists $shortnames{$key}){
	my $orig = $shortnames{$key};
	$shortnames{$key} = $short if(length $orig > length $short);
    } else {
	$shortnames{$key} = $short;
    }
    $lenLong = length $short if(length $short > $lenLong);

    foreach $etype (keys %refEnergies){
	$referenceE{$key}{$etype} = $refEnergies{$etype};
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
	sprintf "..  %-30s%1s%7s %5s  %10.5f %10.5f %3i:%-3i %5i %16.10f %4s %5s",
	"$base","$machine","$dt","$effnw","$hfe","$vare",
	$numci,$numbf,$numjw,$eavg,
	$numerrors,$numwarnings;

    if($wallclock ne ""){
	#the calculation completed, and some extra data is available
	if(abs($corLength) < 1e-10){
	    #The old format of output printed the Sample variance directly
	    $sampleVar = (split/\s+/,`grep "Sample variance" $base.out`)[3];
	    $corLength = $sampleVarCorLen / $sampleVar;
	} else {
	    $sampleVar       = $sampleVarCorLen / $corLength;
	}

	#This is similiar to the Kappa from the 2007 Dolg ECP paper.
	#Lower is better. Sample clock is in microseconds.
	my $wfEfficiency = $dt * $sampleVar * $corLength * $sampleclock * 10.0;

	$summary{$key} .= 
	    sprintf " %10.3e %10.2f %10.2f %10s %10s %10s %15.5f\n",
	    $sampleVar,
	    $corLength,
	    $wfEfficiency,
	    $wallclock,$totalclock,
	    $iteration,($effdt*$iteration);
    } else {
	$summary{$key} .= 
	    sprintf " %10s %10s %10s %10s %10s %10s %15.5f\n",
	    "", "", "", "", "",
	    $iteration,($effdt*$iteration);
    }

    #if we are in enhanced text mode, we need to double escape the "_"
    #$base =~ s/_/\\\\_/g;
    printf DATFILE "#%19s %20s %20s %40s\n", "dt=$dt","$base","E=$eavg","$key";
    print DATFILE "$fordatfile\n\n";
} 
close DATFILE;

chomp($headerLine);
if($extraTag ne ""){
    printf "%30s $headerLine  %15s\n$lastlines"," ",$extraTag;
} else {
    printf "%30s $headerLine \n$lastlines"," ";
}

foreach $key (sort byenergy keys %dt_ave_results)
{
    if(!exists $label{$key})
    {
	$label{$key} = sprintf "%2i", (scalar keys %label) + 1;
    }
}
printf "ID  %-30s %7s %5s  %10s %10s %7s %5s %16s %4s %5s %10s %10s %10s %10s %10s %10s %15s\n",
    "File Name","dt","nw","HF E",
    "Var E","CI:BF",
    "NumJW","Avg E","Err","Warn",
    "Variance",
    "Corr Len","WF Eff","Wall","Total","Iter","effdt*iters";

foreach $sum (sort bydt keys %summary)
{
    $summary{$sum} =~ s/../$label{$sum}/;
    print "$summary{$sum}";
}

if($num_results > 0){
    $ave_result /= $num_results;
    #print "Average result = $ave_result\n";
    $labelLen = $lenLong;
    $labelLen = length "Label" if(length "Label" > $labelLen);
    printf "%5s %*s %10s %20s %5s %7s   %-25s %5s %20s %20s %10s\n",
    "ID",$labelLen,"Label",
    "dt","Ref. Energy","Num","CI:BF","NumJW","NumW","Average","Corr. E.","Weight";
    my %qref;
    my %href;
    my $dtref = 0;
    my $cure = "";
    foreach $key (sort byenergy keys %dt_ave_results)
    {
	my @keydata = split/&/,$key;
	$dt_ave_results{$key} /= $dt_num_results{$key};

	if($dt_nme_results{$key} > 0){
	    $dt_err_results{$key} = sqrt($dt_err_results{$key}/$dt_nme_results{$key});
	} else {

	}

	printf "%5i %*s %10s %20s %5i %3i:%-3i   %-25s %5i %20s %20.10f %10.5f\n",
	$label{$key},
	$labelLen,
	$shortnames{$key},
	"$keydata[1]", "$keydata[0]",
	$dt_num{$key},
	$keydata[5],
	$keydata[2],
	$keydata[3],
	$keydata[4],
	getEnergyWError($dt_ave_results{$key},$dt_err_results{$key}),
	($keydata[0]-$dt_ave_results{$key}),
	$dt_num_results{$key};
    }

    print "\n\n";
    #matrix output
    #the data is sorted according to dt first
    #we calculate the difference for for all results available
    #but we don't compare calculations if dt and energy are different
    print "DMC comparisons:\n";
    my %comparisons;

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
	    my $rOrb = $rowdata[6];
	    my $cOrb = $coldata[6];

	    #the results are not comparable if either is zero
	    next if($rMult == 0 || $cMult == 0);
	    #otherwise we'll get two of every comparison
	    next if($r < $c);
	    #This eliminates a lot of the meaningless comparisons
	    my $orbsMatch = ($rMult * $rOrb == $cMult * $cOrb);
	    next if($orbFilter && !$orbsMatch);

	    #So that we're comparing the difference
	    $cMult *= -1;

	    my $diff = $r*$rMult + $c*$cMult;
	    my $stdR = abs($dt_err_results{$row}*$rMult);
	    my $stdC = abs($dt_err_results{$col}*$cMult);
	    my $diffe = sqrt($stdR*$stdR + $stdC*$stdC);
	    
	    $diff *=  $units;
	    $diffe *= $units;

	    my $comparison = "";
	    my $rStr = getEnergyWError($dt_ave_results{$row},$dt_err_results{$row});
	    my $cStr = getEnergyWError($dt_ave_results{$col},$dt_err_results{$col});

	    if($sumResults == 1){
		$comparison .= sprintf "%3i %3i) %6s",$label{$row},$label{$col},$rowdata[1];
		my $rM = $rMult;
		my $cM = $cMult;
		$rM = " " if($rMult == 1);
		$cM = "- " if($cMult == -1);

		$compType = sprintf " ${rM} %*s +${cM} %*s ",
		$lenLong,$shortnames{$row},
		$lenLong,$shortnames{$col};
		$compType =~ s/\+\-/\-/g;
		$comparison .= sprintf "%s =", $compType;
	    } else {
		$comparison .= sprintf "%3i) %*s %15s %6s %3s:%2s:%-3s %5s | ",
		$label{$row},$lenLong,$shortnames{$row},
		$rStr,$rowdata[1],$rowdata[5],$rOrb,$rowdata[2],$rowJW;
		$comparison .= sprintf "%3i) %*s %15s %6s %3s:%2s:%-3s %5s | ",
		$label{$col},$lenLong,$shortnames{$col},
		$cStr,$coldata[1],$coldata[5],$cOrb,$coldata[2],$colJW;
		$compType = "${rMult}A+${cMult}B";
		$compType =~ s/1//g;
		$compType =~ s/\+\-/\-/g;
		$comparison .= sprintf "%6s=", $compType;
	    }

	    if($orbsMatch){
		$comparison .= " ";
	    } else {
		$comparison .= "*";
	    }

	    $dStr = getEnergyWError($diff,$diffe);
	    if(0){
		$comparison .= sprintf " %9.5f",$diff;
		$comparison .= sprintf " +/- %-9.5f $unitsL", $diffe;
	    } else {
		$comparison .= sprintf " %10s $unitsL",$dStr;
	    }

	    if($compareE){
		foreach $etype (reverse sort keys %{$referenceE{$row}}){
		    $eRow = $referenceE{$row}{$etype} * $rMult;
		    $eCol = $referenceE{$col}{$etype} * $cMult;
		    my $temp = ($eRow + $eCol)*$units;
		    
		    next if(!exists $referenceE{$col}{$etype} || abs($temp) < 1e-10);

		    if($sumResults == 1){
			#$comparison .= sprintf " %*s = %9.5f %s\n",(2*$lenLong+22),"",$temp,$etype;
			$comparison .= sprintf " %s = %9.5f",$etype,$temp;
		    } else {
			$comparison .= sprintf "\n     %*s %15.10f %23s |      %*s %15.10f %23s | %7s  %9.5f %s",
			$lenLong,"",$referenceE{$row}{$etype},"",$lenLong,"",$referenceE{$col}{$etype},"","",$temp,$etype;
		    }
		}
	    }
	    $comparison .= sprintf "\n";

	    if($latexHelp){
		$tempStr = sprintf "%5s & %20s & %15s & %20s & %15s & %15s \\\\",
		$rowdata[1],$shortnames{$row},$rStr,$shortnames{$col},$cStr,$dStr;
		$tempStr =~ s/\./\&/g;
		$comparison = "$tempStr\n";
	    }

	    $comparisons{$comparison} = $diff;
	}
    }

    foreach $key (sort {$comparisons{$a} <=> $comparisons{$b}} keys %comparisons){
	print "$key";
    }
    print "\n\n";
}
