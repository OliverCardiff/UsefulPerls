use strict;
use warnings;

foreach(@ARGV)
{
	my @sps = split(/\./);
	my $stub = $sps[0];

	open FIN, $_ or die $!;
	my $diff = 0;
	my $maxDiff = 0;
	my $comp = 0;
	my $old1 = 0;
	my $old2 = 0;
	my $maxlen = 0;
	my $len = 0;
	while(<FIN>)
	{
		if($. > 1)
		{
			chomp;
			my @sps1 = split(/\t/);
			
			if((scalar @sps1) > 1)
			{
				if($sps1[0] ne "NA")
				{
					$diff = $sps1[0] - $sps1[1];
					if($comp == 1)
					{
						$len = $sps1[0] - $old1;

						if($len > $maxlen)
						{
							$maxlen = $len;
							$maxDiff = $diff;
						}
					}
					
					$comp = 1;
					$old1 = $sps1[0];
					$old2 = $sps1[1];
				}
				else
				{
					$comp = 0;
				}
			}
		}
	}

	close FIN;

	print $stub . "\t" . $maxDiff . "\t" . $maxlen . "\n";
}
