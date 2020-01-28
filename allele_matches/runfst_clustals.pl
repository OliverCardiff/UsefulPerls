use strict;
use warnings;

foreach(@ARGV)
{
	my $fn = $_;
	my @sps = split(/\./);
	
	system("clustalw2fasta $fn" . " > " . $sps[0] . ".fa");
	system("sed -i 's/-/N/g' " . $sps[0] . ".fa");
}
