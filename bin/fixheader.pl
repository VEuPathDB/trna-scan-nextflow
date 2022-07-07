#!/usr/bin/perl

use strict;
use warnings;

my $file = 'fileWithHeader.txt';
open my $info, $file or die "Could not open $file: $!";
open(OUT, ">datafile");

print OUT "Sequence\t\t\ttRNA\tBounds\ttRNA\tAnti\tIntron\tBounds\tInf\nName\t\t\ttRNA#\tBegin\tEnd\tType\tCodon\tBegin\tEnd\tScore\tNote\n--------\t\t------\t-----\t------\t----\t-----\t-----\t----\t------\t------\n"; 
while( my $line = <$info>)  {
    if($line !~ /^Sequence/ and $line !~ /^Name/ and $line !~ /^--------/) {
	print OUT $line;
    }
}
close OUT;
close $info;

