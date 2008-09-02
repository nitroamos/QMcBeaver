#!/usr/bin/perl
#assume utilities.pl is in the same directory as summary.pl
my $path = `dirname $0`;
chomp($path);
require "$path/utilities.pl";

my $printFunc = 1;
my $useScaled = 0;
my $multiPlot = 1;
my $showOpt   = 0;
my $makeGraph = 0;
#put the jastrow in the exponential
my $useExp    = 1;
#square the whole thing (so that the y axis
#can be interpreted as a percentage)
my $useSqr    = 0;

my $date = `date`;
chomp $date;

while($#ARGV >= 0 && $ARGV[0] =~ /^-/){
    $type = shift(@ARGV);
    $param = "";

    if($type eq "-o"){
	$showOpt = ! $showOpt;
	print "Using showOpt = $showOpt\n";
    } elsif($type eq "-p"){
	$makeGraph = 1;
	print "Using makeGraph = $makeGraph\n";
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
    } else {
	print "Unrecognized option: $type\n";
	die;
    }
}

my @files = sort @ARGV;
if($#files < 0){
    push(@files,".");
}

getFileList(".out",\@files);

%jastrows;
%plotters;
my %optEnergies;
my $base = "";
my $numjw = "";
my $numjwID = 0;
my $numbf = 0;
my $numci = 1;
my $refE  = 0;
my $step  = 1;
my $lastS = "";
my $short = "";
for(my $index=0; $index<=$#files; $index++){
    $lastS = $short;
    $base = substr($files[$index],0,-4);
    $short = `basename $base`;
    chomp $short;
    $short =~ s/_[\d]+$//g;

    #print "base = $base\n";
    open (CKMFFILE, "$base.ckmf");
    my $isd = "";
    while(<CKMFFILE>){	
	if($_ =~ m/^\s*run_type\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $isd = $line[1];
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
	if($_ =~ m/^\s*energy\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $refE = $line[1];
	}
	if($_ =~ m/&geometry$/){
	    last;
	}
    }

    $numjw = "";
    $numjwID = 0;

    while(<CKMFFILE> !~ /Jastrow/){}
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

    open (FILE, "$files[$index]");
    my $line = <FILE>;
    my $name = "";
    my $L = 0;
    my $best;
    if($showOpt == 1){
	$best = $step;
    }

    if($lastS ne $short){
	$step = 1;
	if($showOpt == 1){
	    $best = $step;
	}
    }
    while(!eof FILE)
    {
	my $func = "";
	if(($line =~ /Eup/ || $line =~ /Edn/) && $line !~ /parameters/){
	    $name = $line;
	    chomp($name);
	    if($line =~ /Nuclear/ &&
	       ($line =~ /EupE/ || $line =~ /EdnE/)){
		#it's a 3 body jastrow
		while($line !~ /x/){
		    $line = <FILE>;
		}
		chomp;
		my @line = split/[\s]+/,$line;
		$L = $line[4];
	    } else {
		#it's a 2 body jastrow	    
		while($line !~ /x/){
		    $line = <FILE>;
		}
		chomp;
		my @line = split/[\s]+/,$line;
		$L = $line[4];
		$func = <FILE>;
		chomp($func);
	    }

	    $name =~ s/[:()]//g;
	    $name =~ s/Nuclear//;
	    $name =~ s/EupEup/aa/g;
	    $name =~ s/EupEdn/ab/g;
	    $name =~ s/Eup/u/g;
	    
	    $jastrows{"$name&$best&$refE&$numci,$numbf&$numjw&$short"} = "$step&$L&$func&$base";
	}

	if($line =~ /full step/){
	    $step += 2;
	    if($showOpt == 1){
		$best = $step;
	    }
	}
	if($line =~ /Objective Value/ && $line !~ /params/){
	    $optEnergies{"$short&$best"} = (split/ +/,$line)[6];
	}
	$line = <FILE>;	
    }

    close(FILE);
}

printf "%4s %4s %11s %11s %10s %8s %10s %8s   %-30s   %5s   %-s\n",
    "Type","Iter","RefE","VMC E","Corr E","% diff","L (bohr)","% diff","Jastrow","NumBF","File Name";

my $lastL = 0;
my $lastE = 0;
my $lastN = "";
foreach $key (sort a2n3 keys %jastrows)
{
    #$jastrows{"$name&$best&$refE&$numci,$numbf&$numjw&$base"} = "$step&$L&$func";
    ($jName,$best,$refE,$dType,$jType,$short) = split/&/,$key;
    ($step,$L,$func,$base) = split/&/,$jastrows{$key};

    $nrg = $optEnergies{"$short&$best"};

    $corrE = ($refE-$nrg)*627.5095;

    printf "%-4s %4i %11.6f %11.6f", $jName, $step, $refE, $nrg, $corrE;
    if($jName ne $lastN){
	$lastL = $L;
	$lastE = $corrE;
	$lastN = "$jName";
	printf " %10.3f %8s", $corrE, " ";
	printf " %10.5f %8s", $L, " ";
    } else {
	$diffE = $corrE - $lastE;
	if(abs($corrE) > 1e4){
	    printf " %10.1e ", $corrE;
	} else {
	    printf " %10.3f ", $corrE;
	}
	if(abs($diffE) > 1e4){
	    printf "%8.1e", $diffE;
	} else {
	    printf "%8.2f", $diffE;
	}
	printf " %10.5f %8.2f", $L, 100.0*($L - $lastL)/$L;
    }
    $lastL = $L;
    #$lastE = $corrE;

    printf "   %-30s %7s   %-s\n", $jType, $dType, $base;
    
    if(!($func eq "")){
	$plotters{$jName} .= "$jName&$dType&$L&$jType&$func&$nrg&$short#";
    }    
}

die if($makeGraph == 0);

my $gnuplot = "/ul/amosa/bin/gnuplot";
$base =~ s/_[\d]+$//g if(!$showOpt);
my $modbase = $base;
$modbase =~ s/_/\\\\_/g;
my $printedHeader = 0;
foreach $key (reverse sort keys %plotters)
{
    my $filename;
    if($multiPlot){
	$file_name = "jastrows";
    } else {
	$file_name = "$key";
	$printedHeader = 0;
    }
    if($showOpt){
	$file_name .= "_${base}_plot.pdf";
    } else {
	$file_name .= "_plot.pdf";
    }

    if($showOpt){
	$caption .= ", $modbase";
    } else {
	$caption .= ", key L; CI,BF; JW; ID";
    }
    my $xlabel = "r_{ij} (Bohr)";
    my $ylabel = "u_{ij}";

    if($useExp){
	$ylabel = "Exp[$ylabel]";
	if($useSqr){
	    $ylabel = "|$ylabel|^2";
	}
    }

    print "Adding graph of $key to: $file_name\n";

    if(!$printedHeader){

	#`/bin/rm -f $file_name`;
	open(GNUPLOT, "|$gnuplot");
	#open(GNUPLOT, ">gnuplot.gnu");
	print GNUPLOT <<gnuplot_Commands_Done;
#fonts with extensions "ttf" and "dfont" will work
#here is a list of available fonts: Chalkboard Helvetica Times
#Courier Monaco LucidaGrande
#set term gif crop enhanced font 'Monaco' 8
    
#fonts on hive:
#set term gif crop enhanced font 'VeraMono' 8
#set term svg dynamic enhanced font "VeraMono,8"
	
#fonts built into PDFLib Lite:
#Courier, Courier-Bold, Courier-Oblique, Courier-BoldOblique,
#Helvetica, Helvetica-Bold, Helvetica-Oblique, Helvetica-BoldOblique, 
#Times-Roman, Times-Bold, Times-Italic, Times-BoldItalic, Symbol, ZapfDingbats
set term pdf color enhanced font "Courier-Bold,14" linewidth 5 size 17.5,10

set output "$file_name"

#set term png
#set terminal png medium
set size 0.9,1

unset colorbox
show style line
#set logscale y 2
set grid ytics
set mytics
set tics scale 1.5, 0.75

set nokey
set key outside below box Left
set key reverse
#set key noenhanced
set xlabel "$xlabel"
set ylabel "$ylabel"

#set yrange[$y_min:$y_max]
gnuplot_Commands_Done

	if($multiPlot){
	    $numPlots = scalar keys %plotters;
	    $numR = 2;
	    $numC = 2;
	    $numC = 3 if($numPlots > 4);
	    $numR = 3 if($numPlots > 6);
	    die "Too many plots: $numPlots" if($numPlots > 9);
	    print GNUPLOT "set multiplot layout $numR,$numC\n";
	}
    }

    my $caption = "$key";    
    $caption =~ s/Eup/E_{up} /g;
    $caption =~ s/Edn/E_{dn} /g;
    $caption =~ s/Nuclear([\w]+)/$1/g;
    $caption = "$caption Jastrow Functions";
    if($printedHeader){
	print GNUPLOT "set title \"$caption\"\n";
    } else {
	print GNUPLOT "set title \"$caption\\n{/=8${date}}\"\n";
    }
    $printedHeader = 1;

    #$jastrows{"$name&$best&$refE&$numci,$numbf&$numjw&$base"} = "$step&$L&$func&$energy"; 
    my @plots = split/\#/,$plotters{$key};

    my $xmax = 0;
    my $longestJW = 0;

    for(my $i=0; $i<=$#plots; $i++){
	if($useScaled){
	    $xmax = 1;
	} else {
	    my $new = (split/&/,$plots[$i])[2];
	    if( $new > $xmax){
		$xmax = $new;
	    }
	}
    }

    print GNUPLOT "plot [0:$xmax]";
    for(my $i=0; $i<=$#plots; $i++){
	($jName,$dType,$max,$jw,$func,$optE,$example) = split/&/,$plots[$i];
	$jw =~ s/18,//g;
	$jw =~ s/18//g;

	my $title;

	if($showOpt){
	    $title = sprintf "%2i %8.4f", $i, $optE;
	} else {
	    if($max >= 10.0){
		$title = sprintf "%-4.1f;",$max;
	    } else {
		$title = sprintf "%-4.2f;",$max;
	    }
	    #$title = "";
	    $title .= sprintf " %3s; %s; %s", $dType, $jw, $example;
	}
	$title =~ s/_/\\\\_/g;
	#change the font size of the key
	#$title = "{/=10$title}";

	$func =~ s/\^/**/g;
	$func =~ s/ +//g;
	
	#a polynomial might not be completed. it might end with a +)
	$func =~ s/\+\)/\)/g;
	
	#add in the implicit multiplications
	#this line confuses emacs' indentation algorithm...
	$func =~ s/([\d])([x(])/$1*$2/g;
	
	if($useScaled){
	    $max = 1;
	} else {
	    $func =~ s/x/(x\/$max)/g;
	}
	
	if($useExp){
	    $func = "exp($func)";
	    if($useSqr){
		$func = "(${func})**2";
	    }
	}
	
	#my $linetype = int($i / 4) + 1;
	#print "line number $i has type $linetype\n";
	$func = "x > $max ? 1/0 : $func";
	#$func = "x";
	print GNUPLOT " $func title \"$title\""; 

	#print GNUPLOT " [0:$kd[2]] $func title \"$kd[3]\""; 
	if($i != $#plots){
	    print GNUPLOT ",\\"; 
	}
	print GNUPLOT "\n";
    }

    if($multiPlot){
	#In order to only include the key once, we must make sure that the
	#plots are always sorted the same!!!
	print GNUPLOT "set nokey\n";
    } else {
	close(GNUPLOT);
    }

    
#`/bin/rm $_.dat`;
    
#`open $file_name`;
}
if($multiPlot){
    print GNUPLOT "unset multiplot\n";
    close(GNUPLOT);
}
