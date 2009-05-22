#!/usr/bin/perl
#use strict;
my $path = `dirname $0`;
chomp($path);
require "$path/utilities.pl";

my $orbFilter  = 1;
my $calcDiff = 1;
my $useAvg   = 1;
my $withErr  = 0;
my $spacef   = 0.3;
my $i_active = 1;

#absolute energies (=0) or relative (=1) to each other?
my $shift = 1;

#should the x axis be iteration (=0), samples (=1) or time (=2)?
my $xtype = 0;

#add lines with these values:
my @exact_titles;
my @exact;

my $every    = 15;
if($withErr){
    #error lines can be very messy, so decrease the
    #freqency of points
    $every    = 100;
}

my $units = 627.50960803;
my $unitsL = "kcal/mol";

while($#ARGV >= 0 && $ARGV[0] =~ /^-/){
    $type = shift(@ARGV);
    $param = "";

    if($type eq "-s"){
	$calcDiff = ($calcDiff+1)%2;
	print "Using calcDiff = $calcDiff\n";
    } elsif($type eq "-a"){
	$useAvg = ($useAvg+1)%2;
	print "Using useAvg = $useAvg\n";
    } elsif($type eq "-o"){
	$orbFilter = ($orbFilter+1)%2;
	print "Using orbFilter = $orbFilter\n";
    } elsif($type eq "-i"){
	$i_active = ($i_active+1)%2;
	print "Using interactive = $i_active\n";
    } elsif($type eq "-t"){
	$param = shift(@ARGV);
	$xtype = $param;
	print "Using xtype = $xtype\n";
    } elsif($type eq "-f" || $type eq "-space"){
	$param = shift(@ARGV);
	$spacef = $param;	
	print "Using spacef = $spacef\n";
    } elsif($type eq "-e" || $type eq "-every" ){
	$param = shift(@ARGV);
	$every = $param;
	$withErr = 1;
	print "Using every = $every\n";
    } elsif($type eq "-err" || $type eq "-error"){
	$withErr = ($withErr+1)%2;
	print "Using withErr= $withErr\n";
    } elsif($type eq "-exp"){
	$title = shift(@ARGV);
	$nrg   = shift(@ARGV);
	push(@exact_titles,$title);
	push(@exact,$nrg);
	print "Adding line called $title at $nrg\n";
    } elsif($type eq "-x"){
	$param = shift(@ARGV);
	if($param eq "ch2"){
	    push(@exact_titles,"exp");
	    push(@exact,-9.353);	    
	} elsif($param == 1) {
	    push(@exact_titles,"exp");
	    push(@exact,-21.5539);	    
	    push(@exact_titles,"ccsdt");
	    push(@exact,-22.5373);
	} else {
	    print "Unrecognized energy choice: $param\n";
	}
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
	print "Using $unitsL energy units, conversion = $units\n";
    } else {
	print "Unrecognized option: $type\n";
	die;
    }
}
if($#ARGV >= 0){
    print "Unrecognized options: @ARGV\n";
}
my $d = qx! date +%F.%H-%M-%S !;
chomp($d);
my $date = `date`;
chomp $date;

#my $gnuplot = "/usr/local/bin/gnuplot";
my $gnuplot = "/ul/amosa/bin/gnuplot";

#you might need to add this command to your .cshrc
# `setenv GDFONTPATH /Library/Fonts:/System/Library/Fonts`;
# `setenv GDFONTPATH /usr/share/fonts/bitstream-vera/`;

my @gnutype = split/ +/,`$gnuplot -V`;
if($gnutype[1] < 4.3)
{
    #we need the extra features that version 4.3 has
    print "GNUPLOT version = $gnutype[1] is incompatible for executable $gnuplot\n";
    die;
}

sub deleteData
{
    foreach $test (@_){
	#print "Deleting $test\n";
    }

    open(DATA,"plotfile.dat");

    my $newdata = "";
    my $match = 0;
    my $num   = 0;
    my $line = <DATA>;
    while($line)
    {
	if($num < $#_ + 1){
	    foreach $test (@_){
		if($line =~ /$test/){
		    chomp($line);
		    $match = 1;
		}
	    }
	}

	if($match == 1){
	    $line = <DATA>;
	    while($line =~ /\d/)
	    {
		$line = <DATA>;
	    }
	    $line = <DATA>;
	    $match = 0;
	    $num += 1;
	} else {
	    $newdata .= "$line";
	}
	
	$line = <DATA>;
    }

    close(DATA);

    open(NEWDATA,">new_plotfile.dat");
    print NEWDATA "$newdata";
    close(NEWDATA);

    `mv new_plotfile.dat plotfile.dat`;
}

sub operateTwo
{
    my $newdata = "";

    #
    # Operate on two streams:
    # final = $fconst * $fkey + $sconst * $skey
    my ($fconst,$fkey,$sconst,$skey) = @_;
    printf "%10.5f * (%-60s) + %10.5f * (%-60s)\n",$fconst, $fkey, $sconst, $skey;

    open(DATA,"plotfile.dat");

    my $line = <DATA>;
    $line = <DATA> while($line !~ /$fkey$/ && $line !~ /$skey$/);
    if($line =~ /$skey/ && $fkey ne $skey){
	#we found the second key first, so swap
	$temp = $fkey;
	$fkey = $skey;
	$skey = $temp;

	$temp = $fconst;
	$fconst = $sconst;
	$sconst = $temp;
    }

    my @ftitle = split/[ =]+/,$line;
    $line = <DATA>;
    my @first_data;
    while($line =~ /\d/)
    {
	push(@first_data,$line);
	$line = <DATA>;
    }
    #this is the second blank line
    $line = <DATA>;
    #this is header of the next data
    $line = <DATA>;

    while($line !~ /$skey$/)
    {
	$line = <DATA>;
    }
    my @stitle = split/[ =]+/,$line;
    $line = <DATA>;
    my @second_data;
    while($line =~ /\d/)
    {
	push(@second_data,$line);
	$line = <DATA>;
    }

    if($#first_data < $#second_data){
	#it's easier to add the shorter to the longer
	my @temp = @first_data;
	@first_data = @second_data;
	@second_data = @temp;

	@temp = @ftitle;
	@ftitle = @stitle;
	@stitle = @temp;

	$temp = $fkey;
	$fkey = $skey;
	$skey = $temp;

	$temp = $fconst;
	$fconst = $sconst;
	$sconst = $temp;
    }
    
    $first_max = (split/ +/,$first_data[$#first_data])[1];
    $second_max = (split/ +/,$second_data[$#second_data])[1];
    chomp($first_max);
    chomp($second_max);
    $first_min = (split/ +/,$first_data[0])[1];
    $second_min = (split/ +/,$second_data[0])[1];
    chomp($first_min);
    chomp($second_min);
    
    #print "First max = $first_max Second max = $second_max fmin = $first_min smin = $second_min\n";
    #print "Data first = $#first_data Data second = $#second_data\n";
    #print "last = $first_data[$#first_data]";
    my $s = 0;
    my @sl = split/ +/,$second_data[$s];
    my $si = $sl[1];
    
    ($fe,$fw) = split/:/,$ftitle[5];
    ($se,$sw) = split/:/,$stitle[5];
    
    #if the stream hasn't been a weight yet, then initialize it with 1
    $fw = 1.0 if($fw eq "");
    $sw = 1.0 if($sw eq "");

    my $fbase = `basename $ftitle[3]`;
    chomp($fbase);
    #$fbase =~ s/_[\d]+$//g;
    my $sbase = `basename $stitle[3]`;
    chomp($sbase);
    #$sbase =~ s/_[\d]+$//g;

    my $title_new;
    my $new_weight;
    if($fconst * $sconst > 0){
	#we're adding streams
	if(length $fbase < length $sbase){
	    $title_new = "$fbase";
	} else {
	    $title_new = "$sbase";
	}

	#the weight of the product stream will be
	#the sum of the weights from the input streams,
	#each scaled by a constant
	$new_weight = $fw * $fconst + $sw * $sconst;
	$fconst *= $fw/$new_weight;
	$sconst *= $sw/$new_weight;
    } else {
	#we're subtracting streams

	#$title_new = "${fconst}x${fbase}-${sconst}x${sbase}";
	my $ffactor;
	if(abs($fconst+1) < 1e-5){
	    # -1
	    $ffactor = "-";
	} elsif(abs($fconst-1) < 1e-5){
	    # +1
	    $ffactor = "";
	} else {
	    $ffactor = "$fconst";
	}
	my $sfactor;
	if(abs($sconst+1) < 1e-5){
	    # -1
	    $sfactor = "-";
	} elsif(abs($sconst-1) < 1e-5){
	    # +1
	    $sfactor = "";
	} else {
	    $sfactor = "$sconst";
	}

	#normalize the weights now
	$new_weight = $fw+$sw;
	$title_new = "$fbase:${ffactor}A+${sfactor}B";
	$title_new =~ s/A\+\-B/A\-B/;
	$title_new =~ s/-A\+B/B\-A/;
    }

    my $e_new = sprintf "%-.10f", ($fconst * $fe + $sconst * $se);
    #printf " E_New: $fconst * $fe + $sconst * $se = $e_new\n";
    $e_new .= ":$new_weight";

    $newdata .= sprintf "#%19s %20s %20s %40s", "dt=$ftitle[2]","$title_new","E=$e_new","$ftitle[6]=$ftitle[7]";
    #printf         "#%19s %20s %20s %40s\n", "dt=$ftitle[2]","$ftitle[3]","E=$e_new","$ftitle[6]=$ftitle[7]";
    
    for( my $f=0; $f<=$#first_data; $f++){
	@fl = split/ +/,$first_data[$f];
	$fi = $fl[1];
	
	my $new;
	#num samples
	$new = ($fl[1] + $sl[1])/2;
	$newdata .= sprintf "%20s ", $new;

	#energy
	$new = $fconst * $fl[2] + $sconst * $sl[2];
	$newdata .= sprintf "%20.10e ", $new;
	if($f == 0){
	    #printf "Energy: %10.5f * (%-20.10f) + %10.5f * (%-20.10f) = %20.10f\n",$fconst, $fl[2], $sconst, $sl[2], $new;
	}

	#variance
	$new = ($fconst * $fl[3])*($fconst * $fl[3])
             + ($sconst * $sl[3])*($sconst * $sl[3]);
	$newdata .= sprintf "%20.10e ", sqrt($new);
	
	#num samples
	$new = ($fl[4] + $sl[4])/2;
	$newdata .= sprintf "%20s ", $new;

	$newdata .= sprintf "\n";
	#print "Averaging $f:$s  $fi with $si\n";
	
	while($si < $fi && $s <= $#second_data){
	    @sl = split/ +/,$second_data[$s];
	    $si = $sl[1];
	    $s += 1;
	}
    }
    
    close(DATA);

    $newdata .= "\n\n";
    return $newdata;
}

sub averageTwo
{
    #
    # This function will look for two plots that represent equilvalent data and can
    # be averaged.
    #
    
    my @lines = `grep dt plotfile.dat`;
    chomp(@lines);
    foreach(my $fset=0; $fset<$#lines; $fset++){
	@fdata = split/\s+/,$lines[$fset];
	chomp @fdata;
	foreach(my $sset=$fset+1; $sset<=$#lines; $sset++){
	    @sdata = split/\s+/,$lines[$sset];
	    chomp(@sdata);
	    
	    if($fdata[4] eq $sdata[4]){
		print "Average $fdata[2] with $sdata[2]\n"; 
		my $newdata = operateTwo(1.0,$fdata[4],1.0,$sdata[4]);
		deleteData($lines[$fset],$lines[$sset]);

		open(NEWDATA,">>plotfile.dat");
		print NEWDATA "$newdata";
		close(NEWDATA);
		return 1;
	    }
	}
    }
    return 0;
}

sub subtractTwo
{
    my @lines = `grep dt plotfile.dat`;
    my @keys;
    my @titles;
    foreach $line (@lines){
	my @data = split/ +/,$line;
	chomp @data;
	push(@keys,$data[4]);
	push(@titles,$data[2]);
    }
    @keys = sort byenergy @keys;
    
    my %newdata;
    for(my $i=$#keys; $i>=0; $i--){
	$iKey = $keys[$i];
	my @iData = split/&/,$iKey;
	for(my $j=0; $j<$i; $j++){
	    $jKey = $keys[$j];
	    my @jData = split/&/,$jKey;
	    next if(!areComparable($iKey,$jKey));

	    my $temp = 0;
	    ($iMult,$temp,$jMult) = getFormula($iData[2],0,$jData[2],$orbFilter);
	    my $orbsMatch = 0;
	    $orbsMatch = 1 if($iMult * $iData[6] == $jMult * $jData[6]);
	    next if($orbsMatch == 0 && $orbFilter == 1 && $temp == 0);

	    #the results are not comparable if either is zero
	    next if($iMult == 0 || $jMult == 0);

	    #printf "(%2i,%2i) $iMult x %-60s : $jMult x %-60s\n",$i,$j,$iKey,$jKey;
	    #printf "(%2i,%2i) ",$i,$j;
	    $key = "";
	    $key .= ($iMult > $jMult ? $iMult : $jMult);
	    $key .= "x";
	    $key .= ($iMult < $jMult ? $iMult : $jMult);

	    print "Subtracting: $i) $titles[$i] - $j) $titles[$j]\n";
	    $newdata{"$key"} .= operateTwo($iMult,$iKey,-1.0*$jMult,$jKey);
	}
    }

    return 0 if(scalar keys %newdata == 0);

    open(NEWDATA,">new_plotfile.dat");
    foreach $key (reverse sort keys %newdata){
	#assume that one one with the highest numbers in the formula
	#are the ones we want to print
	print NEWDATA "$newdata{$key}";

	#if you want all, commment this line:
	last;
    }
    close(NEWDATA);
    return 0;
}

if($useAvg){
    while(averageTwo()){}
}

if($calcDiff){
    my $once = 0;
    subtractTwo();

    if(-e "new_plotfile.dat"){
	`mv new_plotfile.dat plotfile.dat`;
    }
}

#now it's time to generate gnuplot gifs
my @titles;
my @energies;
my @dt_values;
my @keys;

my $all_dt   = "";
my $all_form = "";

#let's not assume we know what's in the data files
my @lines = `grep dt plotfile.dat`;
chomp(@lines);
foreach $line (@lines)
{
    my @data = split/[= ]+/, $line;
    my ($nrg,$num) = split/:/,$data[5];
    #print "$line\n";
    printf "%-30s: from $num data sets, dt=$data[2], with final energy %20.10e $unitsL\n",$data[3],($nrg*$units);

    if($all_dt eq ""){
	$all_dt = $data[2];
    } elsif($all_dt eq "-1"){
	
    } elsif($data[2] ne $all_dt){
	$all_dt = "-1";
    }

    if($calcDiff){
	my @td = split/:/,$data[3];
	if($all_form eq ""){
	    $all_form = $td[1];
	} elsif($all_form == -1){

	} elsif("$all_form" ne "$td[1]"){
	    $all_form = -1;
	}
    }
}


my $y_min;
my $y_max;
my $y_err;
open (DAT_FILE, "plotfile.dat");
my $line = <DAT_FILE>;
while($line){
    chomp $line;
    my @data = split/[= ]+/, $line;
    
    if($line =~ "dt="){
	push(@energies,"$data[5]"); 
	push(@dt_values,"$data[2]"); 
	push(@keys,"$data[6]"); 

	my @td = split/:/,$data[3];
	my $ti = $td[0];
	if($all_dt == -1 && $data[2] != 0){
	    $ti .= ", dt=$data[2]";	    
	}
	if($all_form == -1 && $data[2] != 0){
	    $ti .= ", $td[1]";	    
	}

	my $key = (split/ +/,$line)[4];
	my @kd = split/&/,$key;

	my $bf  = $kd[2];
	my $jw   = $kd[3];
	
	$jw =~ s/18,//g;
	$jw =~ s/18//g;
		
	$ti .= sprintf " %3s; %s", $bf, $jw;

	push(@titles,"$ti"); 
    }

    $line = <DAT_FILE>;

    #Make sure we have the last line in a series
    if($line !~ /[0-9]/ && "$data[2]" =~ /[0-9]/){
	$y_err = $data[3];
	if($data[2] < $y_min || $y_min == 0){
	    $y_min = $data[2];
	    $y_min -= $y_err if($withErr || $#lines == 0);
	}
	if($data[2] > $y_max || $y_max == 0){
	    $y_max = $data[2];
	    $y_max += $y_err if($withErr || $#lines == 0);
	}
    }
}
close(DAT_FILE);

$y_min *= $units;
$y_max *= $units;
my $intr;
my $reference;
if($calcDiff == 0)
{
    #make a guess for the stoicheometry
    my $ratio = $y_min / $y_max;
    $intr = int($ratio + 0.5);
    #and shift the axis to reflect this
    #we'll have gnuplot shift the plots
    $y_max *= $intr;
    $reference = $y_min;

    if($shift == 1){
	$shift = $y_min;
	$y_max = $y_max-$y_min;
	$y_min = 0;
    } else {
	$shift = 0;
    }
} else {
    # since we use shift as a flag and a parameter, we need to
    # set to zero before plotting
    $shift = 0;
}

my $space = $spacef*($y_min - $y_max);
$y_min += $space;
$y_max -= $space;

my $file_name = "qmc";
my $title_extra = "";
if($all_dt != -1 && $all_dt != 0){
    $title_extra .= ", dt=$all_dt";
    $file_name .= "_$all_dt";
}
if($all_form != -1 && $all_form ne ""){
    $title_extra .= ", $all_form";
}

my $ylabel = "Energy ($unitsL)";
my $xindex;
my $xlabel;
if($xtype == 0){
    $xindex = 4;
    $xlabel = "Num Iterations";
} elsif($xtype == 1){
    $xindex = 1;
    $xlabel = "Num Samples";
} elsif($xtype == 2){
    $xindex = 4;
    $xlabel = "Time (Hartrees^{-1})";
}
$file_name .= sprintf "_%i", ($#{titles}+1);
$gnuplot .= " -geometry 1280x740"; #this is optimized for Amos' laptop...
open(GNUPLOT, "|$gnuplot") or die "Can't open GNUPLOT= $gnuplot\n";
#open(GNUPLOT, ">gnuplot.gnu") or die "Can't open GNUPLOT= $gnuplot\n";

if($i_active){
    print "Plotting graph $file_name with X11\n";
    #print GNUPLOT "set terminal x11 reset persist enhanced font \"Courier-Bold,12\" linewidth 2\n";
    print GNUPLOT "set terminal x11 persist raise enhanced font \"Courier-Bold,12\" title \"$file_name\" dashed linewidth 2\n";
} else {
    $file_name .= sprintf "_$d.pdf", ($#{titles}+1);
    print "Writing graph in: $file_name\n";
    `/bin/rm -f $file_name`;
    print GNUPLOT "set term pdf color enhanced font \"Courier-Bold,12\" linewidth 7 dashed dl 3 size 17.5,10\n";
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
set mouse zoomjump
set size 0.9,1

set nokey
set key outside below box noenhanced Left reverse
set yrange[$y_min:$y_max]
set xrange[0:]
set title "QMC Runs${title_extra}\\n{/=8${date}}"
set xlabel "$xlabel"
set ylabel "$ylabel"
set grid ytics
set mytics
set tics scale 1.5, 0.75

gnuplot_Commands_Done

my $numLC = 11;
my @goodlt;
push(@goodlt,1);
push(@goodlt,3);
push(@goodlt,5);
push(@goodlt,4);
push(@goodlt,6);
push(@goodlt,7);


my $plotline = "plot ";
if($#exact >= 0){
    for(my $i=0; $i<=$#exact; $i++){
	$plotline .= "$exact[$i] title \"$exact_titles[$i]\" with lines,\\\n";
    }
}
for(my $i=0; $i<=$#titles; $i++){
    my $factor = 1;

    if($calcDiff == 0){
	#now we calculate the factor used to indicate stoicheometry
	if( abs($intr * $energies[$i] - $reference/$units) < 0.1 && $intr != 1)
	{
	    $factor *= $intr;
	    $titles[$i] .= " x$factor";
	}
    }

    my $xfactor = 1;
    $xfactor = $dt_values[$i] if($xtype == 2);
    
    my $lt = $goodlt[int($i/$numLC)];    
    my $lc = $i % $numLC;

    $plotline .= " \"plotfile.dat\" index $i every vEvery using (\$$xindex * $xfactor):(\$2*vUnits*$factor-$shift):(\$3*vUnits) lc $lc lt $lt title \"$titles[$i]\"";
    $plotline .= " with yerrorlines";
    #$plotline .= ",\\" if($i != $#titles);
    #$plotline .= "\n";
    $plotline .= "," if($i != $#titles);

}
my $plotline_noerr = $plotline;
$plotline_noerr =~ s/yerrorlines/lines/g;
print GNUPLOT "vEvery = $every\n";
print GNUPLOT "vUnits = $units\n";
if($withErr){
    print GNUPLOT "$plotline\n";
} else {
    print GNUPLOT "$plotline_noerr\n";
}

print GNUPLOT "v=0\n";
print GNUPLOT "bind e 'v=v+1; if(v%2) $plotline; else $plotline_noerr'\n";
print GNUPLOT "k=0\n";
print GNUPLOT "bind k 'k=k+1; if(k%2) set nokey; replot; else set key; replot'\n";
print GNUPLOT "bind '-' 'vEvery=vEvery+5; if(v%2) $plotline; else $plotline_noerr'\n";
print GNUPLOT "bind '=' 'vEvery=vEvery-5; if(vEvery < 1) vEvery = 1; if(v%2) $plotline; else $plotline_noerr'\n";
print GNUPLOT "bind '1' 'vUnits=1; set yrange [$y_min/$units:$y_max/$units]; if(v%2) $plotline; else $plotline_noerr'\n";

print GNUPLOT "pause mouse button2\n";
#print GNUPLOT "pause -1 'Hit return to continue'\n";
#print GNUPLOT "pause -1\n";
close (GNUPLOT);
#`/bin/rm $_.dat`;
#`open $file_name`;
if($i_active == 0){
    `bash -c \"echo Current directory \" | /usr/bin/mutt -s \"[jastrows] $file_name\" -a $file_name nitroamos\@gmail.com`;
    `rm $file_name`;
}
