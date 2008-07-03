#!/usr/bin/perl
# This script will help you examine many ckmf files simultaneously
# If you add a 1 after the executable name, it will reset the reference output files.

my $path = `dirname $0`;
chomp($path);
require "$path/utilities.pl";

die "Usage: $0 <QMcBeaver.x>" if($#ARGV < 0);
my $tol = 1e-7;

my $testmode = 1;
$testmode = 0 if($#ARGV == 1);

my $exe = $ARGV[0];
my $shortexe = `basename $exe`;
chomp($shortexe);

if($testmode == 0){
    print "We'll set the data to $exe\n";
} else {
    print "Comparing saved output with output from $exe, with tolerance $tol\n";
}

my @files;
push(@files,".");
getFileList(".ckmf",\@files);

my $header = getCKMFHeader();
print "$header";

my $base = "";
for(my $index=0; $index<=$#files; $index++){
    $base = substr($files[$index],0,-5);
    next if($files[$index] =~ /.\d\d.ckmf$/ || !(-e "${base}.out"));
    my $oneliner = getCKMFSummary("$files[$index]");
    print "$oneliner";
}

print "\n\n";

for(my $index=0; $index<=$#files; $index++){
    $base = substr($files[$index],0,-5);
    my $basename = `basename $base`;
    my $dirname = `dirname $base`;
    chomp($basename);
    chomp($dirname);
    next if($files[$index] =~ /.\d\d.ckmf$/ || !(-e "${base}.out"));

    my $outname;

    if($testmode == 1){
	$outname = "$dirname/${shortexe}.${basename}.out";
	print "Beginning test for $files[$index]:\n";
    } else {
	$outname = "${base}.out";
    }
    $refCompiler = `grep Compiler ${base}.out`;
    $description = `grep "This will" ${base}.ckmf`;
    my $command = "$exe $files[$index] > $outname";
    print "$description";
    print "$command\n";
    print "Reference compiler $refCompiler";
    `$command`;

    next if($testmode == 0);

    my @origE;
    getEnergies("${base}.out",\@origE);
    my @newE;
    getEnergies($outname,\@newE);
    $newCompiler = `grep Compiler $outname`;
    print "     your compiler $newCompiler";
    my $pass = 1;
    my $avg  = 0;
    my $num  = 0;
    
    if($#origE != $#newE){
	print "Why do they have different numbers of energies ($#origE and $#newE)???\n";
    }

    for($i=0; $i<$#origE; $i++){
	$relerror = abs(($origE[$i] - $newE[$i]) / $origE[$i]);
	$avg += $relerror;
	$num += 1;
	if($relerror > $tol)
	{
	    printf "Original Energy = %-20.10f New Energy = %-20.10f have relative error %20.10e\n",$origE[$i],$newE[$i],$relerror;
	    $pass = 0;
	}
    }
    $avg /= $num if($num > 0);

    if($pass == 0 || $num <= 0 || abs($avg) > $tol){
	printf "Test\n*****----> FAILED <----*****\naverage error %20.10e for $num samples\n", $avg;
    } else {
	printf "Test passed, average error %20.10e for $num samples\n", $avg;
    }
    print "\n";
}
