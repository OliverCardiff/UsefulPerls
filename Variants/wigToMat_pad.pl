use strict;
use warnings;

open FAI, $ARGV[1] or die $!;

my %lenHash;

while(<FAI>)
{
	chomp;
	
	my @sps = split(/\t/);
	
	$lenHash{$sps[0]} = $sps[1];
}

close FAI;


open WIG, $ARGV[0] or die $!;
my $oldID = "";
my $stSwitch = 1;
my $oldBase = 0;
my $oldEnd = 0;

while(<WIG>)
{
	chomp;
	if($. > 1)
	{
		my $subs = substr $_, 0, 1;
		
		if($subs eq "v")
		{
			if($stSwitch == 0 && $oldEnd < ($lenHash{$oldID} - 10))
			{
				my $base = $oldEnd;
				
				while($base < ($lenHash{$oldID} - 10))
				{
					print $oldID . "\t" . "0" . "\n";
					$base += 10;
				}
			}
			#my @sps = split(/_/);
			my @sps = split(/=/);
			
			$oldID = $sps[1];
			$stSwitch = 1;
		} else {
			my @sps = split(/\t/);
			if($stSwitch == 0 && $sps[0] > $oldEnd + 10)
			{	
				my $base = $oldEnd ;
				while($sps[0] > $base + 10)
				{
					$base += 10;
					print $oldID . "\t" . "0" . "\n";
				}
			}
			if($stSwitch == 1 && $sps[0] > 20)
			{
				my $base = $sps[0] % 10;
				$oldBase = $base;
				while($base < $sps[0])
				{
					print $oldID . "\t" . "0" . "\n";
					$base += 10;
				}
			}
			print $oldID . "\t" . $sps[1] . "\n";
			$stSwitch = 0;
			$oldEnd = $sps[0];
		}
		
	}
}

close WIG;
