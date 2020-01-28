use strict;
use warnings;

foreach(@ARGV)
{
	my $fn = $_;
	my @sps = split(/\./);
	
	print "lastz NG1/" . $fn . " NG2/" . $fn . " --output=lzoutNG/" . $sps[0] . ".out --rdotplot=lzoutNG/" . $sps[0] . ".plot\n";
}
