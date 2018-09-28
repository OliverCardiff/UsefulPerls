use strict;
use warnings;

open VCF, $ARGV[0] or die $!;

while(<VCF>)
{
	chomp;
	
	my $subs = substr $_, 0, 1;
	if($subs ne "#")
	{
		my @spsA = split(/\t/);
		my @sps = split(/e/, $spsA[0]);
		
		if($sps[1] > 1000000)
		{
			print $_ . "\n";
		}
		
	}
	else
	{
		print $_ . "\n";
	}
}

close VCF;
