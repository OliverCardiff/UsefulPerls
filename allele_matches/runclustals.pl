use strict;
use warnings;

foreach(@ARGV)
{
	my $fn = $_;
	
	system("clustalw $fn\n");
}
