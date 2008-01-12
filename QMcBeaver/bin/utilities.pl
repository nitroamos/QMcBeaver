#!/usr/bin/perl

#alphabet first, numerical second
sub a1n2 {
    my @adata = split/&/,$a;
    my @bdata = split/&/,$b;
    $bdata[1] <=> $adata[1];
    if($adata[0] eq $bdata[0]){
	if($adata[3] eq $bdata[3]) {
	    $bdata[1] cmp $adata[1];
	} else {
	    $adata[3] <=> $bdata[3];
	}
    } else {
	$bdata[0] cmp $adata[0];
    }
}

sub byenergy {
    my @adata = split/&/,$a;
    my @bdata = split/&/,$b;
    if($adata[0] != $bdata[0]){
	$bdata[0] <=> $adata[0];
    } else {
	$bdata[1] <=> $adata[1];
    }
}

sub bydt {
    my @adata = split/&/,$a;
    my @bdata = split/&/,$b;
    if($adata[1] != $bdata[1]){
	#compare dt
	$bdata[1] <=> $adata[1];
    } else {
	#compare energies
	$bdata[0] <=> $adata[0];
    }
}

sub gcf {
    my ($x, $y) = @_;
    ($x, $y) = ($y, $x % $y) while $y;
    return $x;
}

sub getFormula {
    my ($one, $two) = @_;
    my @od = split/&/,$one;
    my @td = split/&/,$two;
    
    if($od[0] != $td[0]){
	return (0,0) if($od[0] == $td[0] ||
			$od[1] != $td[1] ||
			!($od[3] eq $td[3]));
    } else {
	return (0,0);
    }

    $or = $od[2];
    $tr = $td[2];
    $factor = 100;

    while($factor != 1){
	$factor = gcf($or,$tr);
	#print "gcf($or,$tr) = $factor\n";
	$or /= $factor;
	$tr /= $factor;
    }

    #print "ratio = $ratio : intr = $intr : intdiff = $intdiff\n";
    if($or == int($or) && $tr == int($tr) &&
       $or < 10 && $tr < 10)
    {
	return ($tr, $or);
    } else {
	return (0, 0);
    }
}

# this function will fill in the files array with
# all files that have the extension ext. there are
# a few known directories it will not descend into
sub getFileList
{
    my ($ext, $files) = @_;

    #this will scan through all the subdirectories looking for .out files
    my $clean = 0;
    my $loops = 0;
    while($clean == 0)
    {
	$loops++;
	$clean = 1;
	my @newfiles;
	
	for(my $index=0; $index<=$#$files; $index++)
	{
	    my $cur = ${@$files}[$index];
	    chomp($cur);
	    if(-d $cur && $cur !~ /src$/ && $cur !~ /bin$/ && $cur !~ /include$/ && $cur !~ /hide$/)
	    {
		my @list = `ls $cur`;
		foreach $item (@list)
		{
		    #we have a directory in the list, so we're going to need to loop again
		    $clean = 0;
		    chomp($item);
		    if($cur eq ".")
		    {
			push(@newfiles,"$item");
		    } else {
			push(@newfiles,"$cur/$item");
		    }
		}	    
	    } elsif($cur =~ /$ext$/ && $cur !~ /.step[\d]+./ && $cur !~ /.opt[\d]+./){
		push(@newfiles,$cur);
	    }
	}
	@$files = @newfiles;
		
	if($loops > 8)
	{
	    print "Stopping recursion at $loops.\n";
	}
    }
}
1;
