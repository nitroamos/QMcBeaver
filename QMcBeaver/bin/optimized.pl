#!/usr/bin/perl
#assume utilities.pl is in the same directory as summary.pl
my $path = `dirname $0`;
chomp($path);
require "$path/utilities.pl";

my $publication = 0;
my $printFunc = 1;
my $useScaled = 0;
my $multiPlot = 1;
my $showOpt   = 1;
my $makeGraph = 1;
#my $summary   = 1;
my $i_active  = 1;
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
	$makeGraph = ($makeGraph+1)%2;
	print "Using makeGraph = $makeGraph\n";
    } elsif($type eq "-i"){
	$i_active = ($i_active+1)%2;
	print "Using i_active = $i_active\n";
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
	exit;
    }
}

my @files = sort @ARGV;
if($#files < 0){
    push(@files,".");
}

getFileList(".out",\@files);
$showOpt = 1 if($#files == 0);

my $Cnormal = "\x1b[0m";
my $Chilite = "\x1b[37m";

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
    my @stuff = split/[\s:]+/,getCKMFSummary("$base.ckmf");
    $numci  = $stuff[7];
    $numbf  = $stuff[8];
    $optStr = $stuff[9];
    $refE   = $stuff[10];

    open (CKMFFILE, "$base.ckmf");
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
    my $name = "";
    my $L = 1;
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

    my $iterNRG = 0;
    my $iterSTD = 0;
    my $iterN   = 0;
    my @dat;
    my $line;
    while(!eof FILE)
    {
	$line = <FILE>;	
	@dat = split/\s+/,$line;
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
		@dat = split/\s+/,$line;
		chomp;
		$L = $dat[4];
	    } else {
		#it's a 2 body jastrow	    
		$line = <FILE>;
		my $type = $line;
		chomp;
		if($type =~ /Fixed/){
		    $func = <FILE>;
		    chomp($func);
		    $func .= <FILE>;
		    chomp($func);
		    $L = 10.0;
		} else {
		    $line = <FILE>;
		    $line = <FILE>;
		    chomp($func);
		    @dat = split/\s+/,$line;
		    $L = $dat[4];
		    $func = <FILE>;
		}
		chomp($func);
	    }

	    $name =~ s/[:()]//g;
	    $name =~ s/Nuclear//;
	    $name =~ s/EupEup/Parallel Spin/g;
	    $name =~ s/EupEdn/Opposite Spin/g;
	    $name =~ s/Eup/Electron\-/g;
	    
	    $jastrows{"$name&$best&$refE&$numci,$numbf&$numjw&$short"} = "$step&$L&$func&$base";
	}

	if($line =~ /full step/){
	    $step += 2;
	    if($showOpt == 1){
		$best = $step;
	    }
	}
	
	if($line !~ /[A-Za-df-z]/ && $#dat == 9){
	    $iterNRG = $dat[2];
	    $iterSTD = $dat[3];
	    $iterN   = int($dat[4]/1000 + 0.5);
	    #print "energy line nrg=$iterNRG std=$iterSTD iterN=$iterN\n";
	}

	if($line =~ /Objective Value/ && $line !~ /params/){
	    $optEnergies{"$short&$best"} = "$iterNRG&$iterSTD&$iterN&$optStr";
	}
    }

    if(!defined $optEnergies{"$short&$best"}){
	$optEnergies{"$short&$best"} = "$iterNRG&$iterSTD&$iterN&$optStr"; 
    }

    close(FILE);
}

printf "%15s %4s %11s %11s %7s  %10s %8s %8s %10s %8s   %-30s   %5s    %8s %-s\n",
    "Type","Iter","RefE","VMC E",getOPTHeader(),"Corr E ","std.e.","% diff","L (bohr)","% diff","Jastrow","NumBF","NSmpl(k)","File Name";

my $lastL = 0;
my $lastE = 0;
my $lastN = "";

my @optAvg;
my @optWeight;
my $avgLen = 4;
my $startStep = -1;

foreach $key (sort a2n3 keys %jastrows)
{
    #$jastrows{"$name&$best&$refE&$numci,$numbf&$numjw&$base"} = "$step&$L&$func";
    ($jName,$best,$refE,$dType,$jType,$short) = split/&/,$key;
    ($step,$L,$func,$base) = split/&/,$jastrows{$key};

    ($nrg,$std,$nsamples,$optStr) = split/&/,$optEnergies{"$short&$best"};
    $corrE = ($refE-$nrg)*627.5095;
    $std  *= 627.5095;

    if($base ne $startStep){
	@optAvg    = ();
	@optWeight = ();
    }
    $startStep = $base;

    my $stepVar = 0;
    push(@optAvg,$corrE);
    push(@optWeight,$std);
    shift @optAvg if($#optAvg >= $avgLen);
    shift @optWeight if($#optAvg >= $avgLen);
    my $x  = 0;
    my $x2 = 0;
    my $ws = 0;
    for(my $i=0; $i<=$#optAvg; $i+=1){
	my $val = $optAvg[$i];
	my $w   = $optWeight[$i];
	$ws += $w;
	$x  += $val*$w;
	$x2 += $val*$val*$w;
    }
    $x  /= $ws if(abs($ws) > 0);
    $x2 /= $ws if(abs($ws) > 0);
    $stepVar = $x2 - $x*$x;
    $stepVar = sqrt(abs($stepVar));

    printf "%-15s %4i %11.6f %11.6f %7s", $jName, $step, $refE, $nrg, $optStr;

    my $corrEstr = "";
    if(abs($corrE) > 1e4 || $std == 0){
	$corrEstr = sprintf "  %10.1e", $corrE;
    } else {
	$corrEstr = sprintf "  %-10s", getEnergyWError($corrE,$std);
    }
    if($corrE < 0){
	printf "$Chilite$corrEstr$Cnormal";
    } else {
	printf "$corrEstr";
    }

    printf " %8.2f",$stepVar;
    if($jName ne $lastN){
	$lastL = $L;
	$lastE = $corrE;
	$lastN = "$jName";
	printf " %8s", "";
	printf " %10.5f %8s", $L, " ";
    } else {
	$diffE = $corrE - $lastE;
	if(abs($diffE) > 1e4){
	    printf " %8.1e", $diffE;
	} else {
	    printf " %8.2f", $diffE;
	}
	printf " %10.5f %8.2f", $L, 100.0*($L - $lastL)/$L;
    }
    $lastL = $L;
    #$lastE = $corrE;

    printf "   %-30s %7s    %8g %-s\n", $jType, $dType, $nsamples, $base;
    
    if(!($func eq "")){
	$plotters{$jName} .= "$jName&$dType&$L&$jType&$func&$nrg&$short&$step#";
    }    
}

exit if($makeGraph == 0);

my $gnuplot = "/ul/amosa/bin/gnuplot";
$base =~ s/_[\d]+$//g if(!$showOpt);
my $modbase = $base;
$modbase =~ s/_/\\\\_/g;
my $printedHeader = 0;
my @goodlt;
push(@goodlt,3);# if($publication == 0);
push(@goodlt,1);
push(@goodlt,5);
push(@goodlt,4);
push(@goodlt,6);
push(@goodlt,7);

my $allPlots = "set key outside below box Left reverse;\\\n";
my $lastPlot = "set key outside below box Left reverse;\\\n";

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
	#$file_name .= "_${base}_plot.pdf";
	$file_name .= "_plot.pdf";
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
	if($i_active){
	    $gnuplot .= " -geometry 1280x740"; #this is optimized for Amos' laptop...      
	    open(GNUPLOT, "|$gnuplot");
	    print GNUPLOT "set terminal x11 persist raise enhanced font \"Courier-Bold,12\" title \"$file_name\" dashed linewidth 2\n";
	} else {
	    `/bin/rm -f $file_name`;
	    #open(GNUPLOT, ">gnuplot.gnu");
	    open(GNUPLOT, "|$gnuplot");
	    if($publication == 1){
		print GNUPLOT "set term pdf color enhanced font \"Courier-Bold,16\" linewidth 10 dashed dl 3 size 17.5,10\n";
	    } else {
		print GNUPLOT "set term pdf color enhanced font \"Courier-Bold,14\" linewidth 5 dashed dl 3 size 17.5,10\n";
	    }
	    print GNUPLOT "set output \"$file_name\"\n"; 
	}
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
set size 0.9,1
unset colorbox
show style line
#set logscale y 2
set grid ytics
set mytics
set tics scale 1.5, 0.75
set nokey
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
	    $allPlots .= "set multiplot layout $numR,$numC;\\\n";
	    $lastPlot .= "set multiplot layout $numR,$numC;\\\n";
	}
    }

    my $caption = "$key";    
    $caption =~ s/Eup/E_{up} /g;
    $caption =~ s/Edn/E_{dn} /g;
    $caption =~ s/Nuclear([\w]+)/$1/g;
    $caption = "$caption Jastrow Functions";
    if($printedHeader || $publication == 1){
	$allPlots .= "set title \"$caption\";\\\n";
	$lastPlot .= "set title \"$caption\";\\\n";
    } else {
	$allPlots .= "set title \"$caption\\n{/=8${date}}\";\\\n";
	$lastPlot .= "set title \"$caption\\n{/=8${date}}\";\\\n";
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

    $allPlots .= "plot [0:$xmax]";
    $lastPlot .= "plot [0:$xmax]";
#    for(my $i=0; $i<=$#plots; $i++){
    for(my $i=$#plots; $i>=0; $i -= 1){
	($jName,$dType,$max,$jw,$func,$optE,$example,$step) = split/&/,$plots[$i];
	$jw =~ s/18,//g;
	$jw =~ s/18//g;

	my $title;

	if($showOpt){
	    $title = sprintf "%2i %8.4f", $step, $optE;
	    #$title = $example;
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
	$func =~ s/([\d])([x\(])/$1*$2/g;
	
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
	
	my $lt;
	if($publication == 1){
	    $lt = $goodlt[$i%12];
	} else {
	    $lt = $goodlt[int($i/12)];
	}
	my $lc = ($i+1) % 12;
	#print "line number $i has type lc $lc lt $lt\n";
	$func = "x > $max ? 1/0 : $func";
	#$func = "x";
	$func = " $func lc $lc lt $lt title \"$title\""; 
	$allPlots .= $func;
	
	#print GNUPLOT " [0:$kd[2]] $func title \"$kd[3]\""; 
	if($i == 0){
	    $allPlots .= ";\\"; 
	    $lastPlot .= " $func;\\\n";
	} elsif($i==1) {
	    $allPlots .= ",\\"; 
	    $lastPlot .= " $func,\\\n";
	} else {
	    $allPlots .= ",\\"; 
	}
	$allPlots .= "\n";
    }
    
    if($multiPlot){
	#In order to only include the key once, we must make sure that the
	#plots are always sorted the same!!!
	$allPlots .= "set nokey;\\\n";
	$lastPlot .= "set nokey;\\\n";
    }
    
#`/bin/rm $_.dat`;
    
#`open $file_name`;
}
$allPlots .= "unset multiplot";
$lastPlot .= "unset multiplot";

#print $lastPlot;
#print "bind e 'v=v+1; if(v%2) $lastPlot; else $allPlots;'\n"; 
#die;
print GNUPLOT "v=0\n";
print GNUPLOT "bind l 'v=v+1; if(v%2) $lastPlot; else $allPlots'\n"; 

print GNUPLOT "$allPlots\n";

if($multiPlot){
    #print GNUPLOT "unset multiplot\n";
}
print GNUPLOT "pause mouse button2\n";
close(GNUPLOT);

if($i_active == 0){
    my $email = "nitroamos\@gmail.com";
    print "Check $email...\n";
    `bash -c \"echo Current directory \" | /usr/bin/mutt -s \"[jastrows] $file_name\" -a $file_name $email`;
    `rm $file_name`;
}
