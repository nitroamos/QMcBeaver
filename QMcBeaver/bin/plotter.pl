#!/usr/bin/perl
#use strict;
my $path = `dirname $0`;
chomp($path);
require "$path/utilities.pl";

my $calcDiff = 1;
my $useAvg   = 1;
my $withErr  = 1;

#add lines with these values:
#my @exact_titles    = ("exp" , "ccsdt");
#my @exact           = (-21.5539 , -22.5373);
my @exact_titles    = ("exp");
my @exact           = (-9.353);

my $every    = 1;
if($withErr){
    #error lines can be very messy, so decrease the
    #freqency of points
    $every    = 100;
}
#hartrees (=1) or kcal/mol (=627.50960803)?
my $units = 627.50960803;

#absolute energies (=0) or relative (=1) to each other?
my $shift = 1;

#should the x axis be samples (=1) or iterations (=0)?
my $xaxis_samples = 0;
my $mult_dt       = 1;

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
    $line = <DATA> while($line !~ /$fkey/ && $line !~ /$skey/);

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

    while($line !~ /$skey/)
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
    #print "num first = $#first_data snum = $#second_data\n";
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
    $fbase =~ s/_[\d]+$//g;
    my $sbase = `basename $stitle[3]`;
    chomp($sbase);
    $sbase =~ s/_[\d]+$//g;

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
	$fkey = (split/ +/,$lines[$fset])[4];
	chomp($fkey);
	foreach(my $sset=$fset+1; $sset<=$#lines; $sset++){
	    $skey = (split/ +/,$lines[$sset])[4];
	    chomp($skey);
	    
	    if($fkey eq $skey){
		my $newdata = operateTwo(1.0,$fkey,1.0,$skey);
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
    foreach $line (@lines){
	my $key = (split/ +/,$line)[4];
	chomp($key);
	push(@keys,$key);
    }
    @keys = sort byenergy @keys;

    for(my $i=0; $i<=$#keys; $i++){
	#print "key $i = $keys[$i]\n";
    }
    
    my %newdata;
    for(my $i=$#keys; $i>=0; $i--){
	$iKey = $keys[$i];
	for(my $j=0; $j<$i; $j++){
	    $jKey = $keys[$j];

	    ($iMult,$jMult) = getFormula($iKey,$jKey);

	    #the results are not comparable if either is zero
	    next if($iMult == 0 || $jMult == 0);

	    #printf "(%2i,%2i) $iMult x %-60s : $jMult x %-60s\n",$i,$j,$iKey,$jKey;
	    #printf "(%2i,%2i) ",$i,$j;
	    $key = "";
	    $key .= ($iMult > $jMult ? $iMult : $jMult);
	    $key .= "x";
	    $key .= ($iMult < $jMult ? $iMult : $jMult);

	    $newdata{"$key"} .= operateTwo($iMult,$iKey,-1.0*$jMult,$jKey);
	}
    }

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
    #print "$line\n";
    my @data = split/[= ]+/, $line;

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
	if($data[2] < $y_min || $y_min == 0){
	    $y_min = $data[2];
	    $y_min -= $data[3] if($withErr);
	}
	if($data[2] > $y_max || $y_max == 0){
	    $y_max = $data[2];
	    $y_max += $data[3] if($withErr);
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

my $space = 0.2*($y_min - $y_max);
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

my $ylabel = "Energy";
if($units == 1){
    $ylabel .= " (au)";
} else {
    $ylabel .= " (kcal/mol)";
}

my $xindex = 4;
my $xlabel = "Num Iterations";
if($xaxis_samples == 1){
    $xindex = 1;
    $xlabel = "Num Samples";
} elsif($mult_dt){
    $xlabel = "Time (Hartrees^{-1})";
}

$file_name .= "_$#{titles}_$d.pdf";
#$file_name .= "_$#{titles}.pdf";

print "Writing graph in: $file_name\n";
`/bin/rm -f $file_name`;
open(GNUPLOT, "|$gnuplot");
print GNUPLOT <<gnuplot_Commands_Done;
#fonts with extensions "ttf" and "dfont" will work
#here is a list of available fonts: Chalkboard Helvetica Times
#Courier Monaco LucidaGrande
#set term gif crop enhanced font 'Monaco' 8

#fonts on hive:
#set term gif crop enhanced font 'VeraMono' 8
#set term svg dynamic enhanced font "VeraMono,8"
set term pdf color enhanced font "Courier-Bold,12" linewidth 7 size 12,8
set output "$file_name"

#set term png
#set terminal png medium
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

print GNUPLOT "plot ";
if($#exact >= 0){
    for(my $i=0; $i<=$#exact; $i++){
	print GNUPLOT "$exact[$i] title \"$exact_titles[$i]\" with lines,\\\n";
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
    $xfactor = $dt_values[$i] if($mult_dt);

    my $plotline = " \"plotfile.dat\" index $i every $every using (\$$xindex * $xfactor):(\$2*$units*$factor-$shift):(\$3*$units) title \"$titles[$i]\"";
    if($withErr){
	print GNUPLOT "$plotline with yerrorlines";
    } else {
	print GNUPLOT "$plotline with lines";
    }
    print GNUPLOT ",\\" if($i != $#titles);
    print GNUPLOT "\n";
}

close (GNUPLOT);
#`/bin/rm $_.dat`;

#`open $file_name`;
