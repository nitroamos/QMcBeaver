#!/usr/bin/perl
#assume utilities.pl is in the same directory as summary.pl
my $path = `dirname $0`;
chomp($path);
require "$path/utilities.pl";

my $printFunc = 1;
my $useScaled = 0;
my $useVar    = 0;
my $multiPlot = 1;

#put the jastrow in the exponential
my $useExp    = 1;
#square the whole thing (so that the y axis
#can be interpreted as a percentage)
my $useSqr    = 0;

my $date = `date`;
chomp $date;

my $showOpt;
my @files = @ARGV;
if($#files < 0)
{
    push(@files,".");
} elsif($#files == 0 && $files[0] =~ /out/){
    $showOpt = 1;
    $useVar  = 1;
    print "Showing optimization steps\n";
}

if($showOpt != 1){
    getFileList(".out",\@files);
}

%jastrows;
%plotters;
my @optEnergies;
my $base = "";
my $numjw = "";
my $numjwID = 0;
my $numbf = 0;
my $numci = 1;
for(my $index=0; $index<=$#files; $index++){
    $base = substr($files[$index],0,-4);

    next if($base !~ /tw/);

    #print "base = $base\n";
    open (CKMFFILE, "$base.ckmf");
    my $isd = "";
    while(<CKMFFILE>){	
	if($_ =~ m/^\s*run_type\s*$/){
	    $_ = <CKMFFILE>;
	    chomp;
	    my @line = split/[ ]+/;
	    $isd = $line[1];
	    if($useVar == 0)
	    {
		last if($isd eq "variational");
	    }
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
	if($_ =~ m/&geometry$/){
	    last;
	}
    }

    if($useVar == 0){
	next if($isd eq "variational");
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
    my $step = 0;
    while( ($line !~ /TheMan/ && $line !~ /There/) ||
	   ($showOpt == 1))
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
		#printf "%-30s %-20s %-20s\n",$files[$index],$name,$L;
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
		#printf "%-30s %-20s %-20s\n",$files[$index],$name,$L;
	    }

	    $name =~ s/[:()]//g;
	    my $val = "";
	    if($showOpt){
		$val = "$name&$numci,$numbf&$L&$step&$func";
	    } else {
		$val = "$name&$numci,$numbf&$L&$numjw&$func";
	    }
	    $jastrows{$val} = $base;
	}

	if($line =~ /full step/){
	    $step += 1;
	}
	if($line =~ /Obj/ && $line !~ /params/){
	    $energy = (split/ +/,$line)[6];
	    push(@optEnergies,$energy);
	    print "Energy is $energy, num = $#optEnergies\n";
	}
	$line = <FILE>;
	last if( eof FILE );
    }

    close(FILE);
}

if($showOpt){
    printf "Jastrow = $numjw; NumCI = $numci; NumBF = $numbf\n";
    printf "%-30s %-15s %6s %-30s %5s   %-40s\n",
    "Type","L (bohr)","% diff","Step","NumBF","File Name";
} else {
    printf "%-30s %-15s %6s %-30s %5s   %-40s\n",
    "Type","L (bohr)","% diff","Jastrow","NumBF","File Name";
}

my $lastL = 0;
foreach $key (sort a1n2 keys %jastrows)
{
    my @keydata = split/&/,$key;
    printf "%-30s %-15s", $keydata[0], $keydata[2];

    printf " %6.2f", 100.0*($keydata[2] - $lastL)/$keydata[2];
    $lastL = $keydata[2];

    printf " %-30s %5s   %-40s", $keydata[3], $keydata[1], $jastrows{$key};

    #printf " $keydata[4]";
    printf "\n";
    
    my $example = `basename $jastrows{$key}`;
    chomp($example);

    #remove any run index identifiers
    $example =~ s/_[\d]+$//g;

    if(!($keydata[4] eq "")){
	$plotters{$keydata[0]} .= "$key&$example#";
    }    
}

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

    print "Writing graph in: $file_name\n";

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
	    print GNUPLOT "set multiplot layout 2,2\n";
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
	my @kd = split/&/,$plots[$i];
	my $func = $kd[4];
	my $max = $kd[2];
	my $jw = $kd[3];

	$jw =~ s/18,//g;
	$jw =~ s/18//g;

	my $title;

	if($showOpt){
	    $title = sprintf "%2i %8.4f", $i, $optEnergies[$i];
	} else {
	    if($max >= 10.0){
		$title = sprintf "%-4.1f;",$max;
	    } else {
		$title = sprintf "%-4.2f;",$max;
	    }
	    $title = "";
	    $title .= sprintf " %3s; %s; %s", $kd[1], $jw, $kd[5];
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
	
	my $linetype = int($i / 4) + 1;
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
