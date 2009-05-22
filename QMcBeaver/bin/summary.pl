#!/usr/bin/perl

# Quick start guide is found by running: summary.pl -h
#
#
#
#use strict;
#assume utilities.pl is in the same directory as summary.pl
my $path = `dirname $0`;
chomp($path);
require "$path/utilities.pl";

# First, select the default values for all our parameters.
my $useVar       = 0;
my $dtFilter     = 0;
my $orbFilter    = 1;
my $compareE     = 0;
my $sumResults   = 1;
my $latexHelp    = 0;
my $latexDTcol   = 0;
my $averageTitle = 0;
my $estd_stop    = 0.0;#in $units
#my $extraTag  = "trail_eps2";

my @fileFilters;
my @exclusionFilters;


my $units = 627.50960803;
my $unitsL = "kcal/mol";

#keep only 1 line every $drop lines
#also look at $every in plotter.pl
my $drop = 1;
if($drop != 1){
    print "Keeping only 1 line in $drop\n";
}

#Second, read in user input.
my @files;
while($#ARGV >= 0){
    $type = shift(@ARGV);
    $param = "";

    if($type !~ /^-/){
	#assume for now that it an output file
	push(@files,$type);
    } elsif($type eq "-h"){
	print "Usage:\n";
	print "-h Print this help.\n";
	print "-v Include VMC calculations (currently = $useVar).\n";
	print "-a Average equivalent files (currently = $averageTitle).\n";
	print "-t <param> Only include dt=<param> (or all if 0, currently = $dtFilter).\n";
	print "-f <param> Only include files that match <param>.\n";
	print "-x <param> Exclude files that match <param>.\n";
	print "-u <param> Convert energy units to <param> units. E.g. <param> = ev or kcal\n";
	print "-o Include comparisons between inconsistent orbitals (currently = $orbFilter).\n";
	print "-c Include non-DMC energy comparisons, if available (currently = $compareE).\n";
	print "-e <param> Stop reading calculations when the error goes below <param>, in the selected units (currently = $estd_stop).\n";
	print "-s Summarize output if sumResults=1 (currently = $sumResults).\n";
	print "-l Make a LaTeX table (currently = $latexHelp).\n";
	print "Any option not starting with a '-' will be interpreted as a calculation file/directory.\n";
	print "Directories are recursively scanned, ignoring any directories named \"hide\"\n";
	print "If you don't include any calculation files, then we'll add directory \".\"\n";
	print "So far, you've selected \"@files\".\n";
	exit;
    } elsif($type eq "-v"){
	$useVar = ($useVar+1)%2;
	print "Using useVar = $useVar\n";
    } elsif($type eq "-a"){
	$averageTitle = ($averageTitle+1)%2;
	print "Using averageTitle = $averageTitle\n";
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
	$param = lc($param);
	if($param =~ /kcal/){
	    $units  = 627.50960803;
	    $unitsL = "kcal/mol";
	} elsif($param =~ /ev/) {
	    $units  = 27.211399;
	    $unitsL = "eV"; 
	} elsif($param =~ /kj/) {
	    $units  = 2625.5002;
	    $unitsL = "kJ/mol"; 
	} elsif($param =~ /cm/) {
	    $units  = 219474.63;
	    $unitsL = "cm^-1"; 
	} elsif($param =~ /au/ || $param =~ /hart/) {
	    $units  = 1;
	    $unitsL = "au"; 
	}
	print "Converting energy units: 1.0 $unitsL = $units au\n";
    } elsif($type eq "-o"){
	$orbFilter = ($orbFilter+1)%2;
	if($orbFilter == 1){
	    print "Filtering to only include balanced orbitals\n";
	} else {
	    print "Not filtering results based on orbital usage.\n";
	}
    } elsif($type eq "-c"){
	$compareE = ($compareE+1)%2;
	print "Comparing with reference energies, compareE = $compareE.\n";
    } elsif($type eq "-e"){
	$param = shift(@ARGV);
	$estd_stop = 1*$param;
	#Print the message later, once we're sure $unitsL has been set
    } elsif($type eq "-s"){
	$sumResults = ($sumResults+1)%2;
	print "Summarize report, sumResults = $sumResults.\n";
    } elsif($type eq "-l"){
	$latexHelp = ($latexHelp+1)%2;
	print "LaTex Helper, latexHelp = $latexHelp.\n";
    } else {
	print "Unrecognized option: $type\n";
	die;
    }
}

if($estd_stop > 0.0){
    print "Notice: we will stop reading calculations once they reach an error of $estd_stop $unitsL!\n\n";
}

push(@files,".") if($#files < 0);
#getFileList(".out",\@files);
getFileList(".qmc",\@files);
open (DATFILE, ">plotfile.dat");

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

    if($averageTitle == 1){
	#remove any _\d at the end
	$short = $1 if($short =~ /([\w\d]+)_[\d]+$/);
    }

    my $vare = "";

    my $dt    = 0;
    my $oepi  = 0;
    my $nw    = 0;
    my $effnw = 0;
    my $opt   = -1;
    my $isd   = -1;
    my $hfe   = 0;
    my $numbf = 0;
    my $numci = 0;
    my $numor = 0;
    my $use3  = 0;
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
	if($_ =~ m/^\s*one_e_per_iter\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $oepi = $line[1];
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

    #next if($opt == 1);
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

	next if($estd > 0 && $estd * $units < $estd_stop);
	#this is to avoid processing lines with warnings
	next if( $_ =~ /[=:]/ && $_ !~ /Results/);

	chomp;
	@data = split/[ ]+/;

        #this is the number of data elements per line
	#it can have the letter 'e' or 'E' since scientific notation uses them
	if($#data >= 8 && $_ !~ /[A-DF-Za-df-z]+/ && $more){
	    $counter++;
	    $iteration   = $data[1];
	    $iteration  /= 8 if($oepi == 1);
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
    #$numwarnings /= $num_samples;
    $numwarnings = sprintf "%4i", $numwarnings;
    $numwarnings = " $Chilite$numwarnings$Cnormal" if($numwarnings > 0.5);
    my $numerrors = `grep ERROR $base.out | wc -l`;
    #$numerrors /= $num_samples;
    $numerrors = sprintf "%4i",$numerrors;

    $lastlines .= "$line";
    my $key = "$dt&$numbf&$numjw&$nw&$numci&$numor&$oepi&$short";
    if($vare eq ""){
	#use the value for energy in the key
	$key =  "$hfe&$key";
    } else {
	#use the variational energy from the header in the key
	$key = "$vare&$key";
    }

    my $runage = getFileAge("$base.out",1);
    # updated in the last 15 minutes
    $short = "*$short" if($runage < 900 && $estd * $units > $estd_stop);
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
	sprintf "..  %-30s%1s%7s %5s %16s %4s %5s",
	"$base","$machine","$dt","$effnw",getEnergyWError($eavg,$estd),$numerrors,$numwarnings;

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
printf "ID  %-30s %7s %5s %16s %4s %5s %10s %10s %10s %10s %10s %10s %15s\n",
    "File Name","dt","nw","Avg E","Err","Warn",
    "Variance",
    "Corr Len","WF Eff","Wall","Total","Iter","effdt*iters";

foreach $sum (sort bydt keys %summary)
{
    $summary{$sum} =~ s/../$label{$sum}/;
    print "$summary{$sum}";
}

die if($num_results <= 0);

$ave_result /= $num_results;
#print "Average result = $ave_result\n";
$labelLen = $lenLong;
$labelLen = length "Label" if(length "Label" > $labelLen);
printf "%5s %*s %10s %1s %20s %5s %7s   %-25s %5s %20s %20s %10s\n",
    "ID",$labelLen,"Label",
    "dt","e","Ref. Energy","Num","CI:BF","NumJW","NumW","Average","Corr. E.","Weight";
my %qref;
my %href;
my $dtref = 0;
my $cure = "";
foreach $key (sort byenergy keys %dt_ave_results)
{
    my @keydata = split/&/,$key;

    if($dt_num_results{$key} > 0){
	$dt_ave_results{$key} /= $dt_num_results{$key};
    } else {
	print "Why does $key have $dt_num_results{$key} results?\n";
	die;
    }

    if($dt_nme_results{$key} > 0){
	$dt_err_results{$key} = sqrt($dt_err_results{$key}/$dt_nme_results{$key});
    } else {
	
    }
    
    printf "%5i %*s %10s %1i %20s %5i %3i:%-3i   %-25s %5i %20s %20.10f %10.5f\n",
    $label{$key},
    $labelLen,
    $shortnames{$key},
    "$keydata[1]", $keydata[7], "$keydata[0]",
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
my %comparisons;

#A + B = C + D
foreach $A (sort bydt keys %dt_ave_results)
{
    my @Adata = split/&/,$A;
    foreach $C (sort bydt keys %dt_ave_results)
    {
	my @Cdata = split/&/,$C;	
	next if(!areComparable($A,$C));

	foreach $B (sort bydt keys %dt_ave_results)
	{
	    #next if($A eq $B || $A eq $C || $B eq $C);
	    next if(!areComparable($A,$B));

	    my @Bdata = split/&/,$B;
	    my $a = $dt_ave_results{$A};
	    my $b = $dt_ave_results{$B};
	    my $c = $dt_ave_results{$C};
	    next if($a < $c || $a < $b); #otherwise we'll get two of every comparison
	    #next if($a < $c); #otherwise we'll get two of every comparison
	    
	    my $aOrb = $Adata[6];
	    my $bOrb = $Bdata[6];
	    my $cOrb = $Cdata[6];	    

	    ($aMult,$bMult,$cMult) = getFormula($Adata[2],$Bdata[2],$Cdata[2],$orbFilter);
	    #print "$a $b $c ($aMult,$bMult,$cMult) \n" if($bMult != 0);
	    next if($a < $b && $bMult > 0);	    
	    next if($aMult == 0 || $cMult == 0); #the results are not comparable if either is zero
	    #This eliminates a lot of the meaningless comparisons
	    my $orbsMatch = 0;
	    $orbsMatch = 1 if($aMult * $aOrb == $cMult * $cOrb);
	    next if($orbsMatch == 0 && $orbFilter == 1 && $bMult == 0);
            #So that we're comparing the difference
	    $cMult *= -1;

	    #print "$orbsMatch = ($aMult * $aOrb == $cMult * $cOrb)   ($aMult,$bMult,$cMult)\n";
	    #print "$orbsMatch \n";


	    my $diff = $a*$aMult + $b*$bMult + $c*$cMult;
	    my $stdA = abs($dt_err_results{$A}*$aMult);
	    my $stdB = abs($dt_err_results{$B}*$bMult);
	    my $stdC = abs($dt_err_results{$C}*$cMult);
	    my $diffe = sqrt($stdA*$stdA + $stdB*$stdB + $stdC*$stdC);
	    
	    $diff  *= $units;
	    $diffe *= $units;
	    
	    my $comparison = "";
	    my $aStr = getEnergyWError($dt_ave_results{$A},$dt_err_results{$A});
	    my $bStr = getEnergyWError($dt_ave_results{$B},$dt_err_results{$B});
	    my $cStr = getEnergyWError($dt_ave_results{$C},$dt_err_results{$C});
	    my $diffStr = getEnergyWError($diff,$diffe);

	    #print "($Adata[2],$Bdata[2],$Cdata[2]) => ($aMult,$bMult,$cMult) := $diffStr\n";
	    
	    if($sumResults == 1){
		$comparison .= sprintf "%3i %3i) %6s",$label{$A},$label{$C},$Adata[1];
		my $aM = $aMult;
		my $bM = $bMult;
		my $cM = $cMult;
		$aM = " " if($aMult == 1);
		$bM = " " if($bMult == 1);
		$cM = "- " if($cMult == -1);
		
		my $compType = sprintf " ${aM} %*s ",
		$lenLong,$shortnames{$A};
		
		if($bMult != 0){
		    $compType .= sprintf " +${bM} %*s ",
		    $lenLong,$shortnames{$B};
		}
		$compType .= sprintf " +${cM} %*s ",
		$lenLong,$shortnames{$C};

		$compType =~ s/\+\-/\-/g;
		$comparison .= sprintf "%s =", $compType;
	    } else {
		my $AJW = (split/=/,$Adata[3])[0];
		my $BJW = (split/=/,$Bdata[3])[0];
		my $CJW = (split/=/,$Cdata[3])[0];
		$comparison .= sprintf "%3i) %*s %15s %6s %3s:%2s:%-3s %5s | ",
		$label{$A},$lenLong,$shortnames{$A},
		$aStr,$Adata[1],$Adata[5],$aOrb,$Adata[2],$AJW;

		if($bMult != 0){
		    $comparison .= sprintf "%3i) %*s %15s %6s %3s:%2s:%-3s %5s | ",
		    $label{$B},$lenLong,$shortnames{$B},
		    $bStr,$Bdata[1],$Bdata[5],$bOrb,$Bdata[2],$BJW;

		    $compType = "${aMult}A+${bMult}B+${cMult}C";
		} else {
		    $compType = "${aMult}A+${cMult}B";
		}

		$comparison .= sprintf "%3i) %*s %15s %6s %3s:%2s:%-3s %5s | ",
		$label{$C},$lenLong,$shortnames{$C},
		$cStr,$Cdata[1],$Cdata[5],$cOrb,$Cdata[2],$CJW;

		$compType =~ s/1//g;
		$compType =~ s/\+\-/\-/g;
		$comparison .= sprintf "%6s=", $compType;
	    }
	    
	    if($orbsMatch){
		$comparison .= " ";
	    } else {
		$comparison .= "*";
	    }
	    
	    if(0){
		$comparison .= sprintf " %9.5f",$diff;
		$comparison .= sprintf " +/- %-9.5f $unitsL", $diffe;
	    } else {
		$comparison .= sprintf " %10s $unitsL",$diffStr;
	    }
	    
	    if($compareE){
		foreach $etype (reverse sort keys %{$referenceE{$A}}){
		    $eA = $referenceE{$A}{$etype} * $aMult;
		    $eB = $referenceE{$B}{$etype} * $bMult;
		    $eC = $referenceE{$C}{$etype} * $cMult;
		    my $temp = ($eA + $eB + $eC)*$units;
		    
		    next if(!exists $referenceE{$C}{$etype} || abs($temp) < 1e-10);
		    
		    if($sumResults == 1){
			#$comparison .= sprintf " %*s = %9.5f %s\n",(2*$lenLong+22),"",$temp,$etype;
			$comparison .= sprintf " %s = %9.5f",$etype,$temp;
		    } else {
			$comparison .= sprintf "\n     %*s %15.10f %23s |      %*s %15.10f %23s | %7s  %9.5f %s",
			$lenLong,"",$referenceE{$A}{$etype},"",$lenLong,"",$referenceE{$C}{$etype},"","",$temp,$etype;
		    }
		}
	    }
	    $comparison .= sprintf "\n";
	    
	    if($latexHelp){
		#strip off any of the title after the first underscore
		my $nameA = $1 if($shortnames{$A} =~ /([\dA-Za-z]+)_/);
		my $nameC = $1 if($shortnames{$C} =~ /([\dA-Za-z]+)_/);

		$tempStr = "";
		$tempStr = sprintf "%5s & ", $Adata[1] if($latexDTcol);		
		$tempStr .= sprintf "%20s & %15s & %20s & %15s & %15s \\\\",
		$nameA,$aStr,$nameC,$cStr,$diffStr;
		$tempStr =~ s/\./\&/g;
		$comparison = "$tempStr\n";
	    }
	    
	    $comparisons{$comparison} = $diff if(abs($diff) < 1000);
	}
    }
}

if($latexHelp){
    if($latexDTcol){
	print <<LATEX_HEADER;
\\begin{center}
\\begin{table}[htdp]
\\caption{A \$\\leftarrow\$ B}
\\label{table:gen}
\\begin{tabular}{r\@{.}l r r\@{.}lr r\@{.}lr\@{.}l}
\\hline \\hline
\\multicolumn{2}{c}{\$\\delta t\$} & A &
\\multicolumn{2}{c}{DMC} & B &
\\multicolumn{2}{c}{DMC} &
\\multicolumn{2}{c}{\$\\Delta\$} \\\\
\\multicolumn{2}{c}{au\$^{-1}\$} &  &
\\multicolumn{2}{c}{au} &  &
\\multicolumn{2}{c}{au} &
\\multicolumn{2}{c}{$unitsL} \\\\
\\hline
LATEX_HEADER
} else {
    print <<LATEX_HEADER;
\\begin{center}
\\begin{table}[htdp]
\\caption{A \$\\leftarrow\$ B}
\\label{table:gen}
\\begin{tabular}{r r\@{.}lr r\@{.}lr\@{.}l}
\\hline \\hline
A &
\\multicolumn{2}{c}{DMC} & B &
\\multicolumn{2}{c}{DMC} &
\\multicolumn{2}{c}{\$\\Delta\$} \\\\
  &
\\multicolumn{2}{c}{au} &  &
\\multicolumn{2}{c}{au} &
\\multicolumn{2}{c}{$unitsL} \\\\
\\hline
LATEX_HEADER
}
}
foreach $key (sort {$comparisons{$a} <=> $comparisons{$b}} keys %comparisons){
    print "$key";
}
if($latexHelp){
print <<LATEX_TAIL;
\\hline \\hline
\\end{tabular}
\\end{table}
\\end{center}
LATEX_TAIL
}

print "\n\n";

