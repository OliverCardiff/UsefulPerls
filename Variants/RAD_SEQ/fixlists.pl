use strict;
use warnings;

open FILE, $ARGV[0] or die $!;

while(<FILE>)
{
	chomp;
	
	my $subs = substr $_, 0, 1;

	if($subs ne "#")
	{	

		my @sps = split(/\t/);
		my @sps2 = split(/\|/, $sps[0]);
		$sps[0] = $sps2[0];
		print join("\t", @sps) . "\n";
	}
}

close FILE;
