#!/usr/bin/perl
#this script will help you examine many ckmf files simultaneously

my $path = `dirname $0`;
chomp($path);
require "$path/utilities.pl";

my $includeRestarts = 0;

my @files = @ARGV;
if($#files < 0)
{
    push(@files,".");
}
getFileList(".ckmf",\@files);

my $header = getCKMFHeader();
print "$header";

my $base = "";
for(my $index=0; $index<=$#files; $index++){
    $base = substr($files[$index],0,-5);

    if(!$includeRestarts){
	next if($files[$index] =~ /.\d\d.ckmf$/ && !(-e "${base}.out"));
    }

    my $oneliner = getCKMFSummary("$files[$index]");
    print "$oneliner";
}
