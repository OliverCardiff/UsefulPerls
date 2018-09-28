use strict;
use warnings;

open BUSCO, $ARGV[0] or die $!;

my %hitHash;

my $oldID;

while(<BUSCO>)
{
	chomp;
	my @sps = split(/\t/);
	
	if($sps[1] ne "Missing")
	{
		push(@{$hitHash{$sps[0]}}, $sps[2]);
	}
	
}

close BUSCO;

my %idxHash;

open IDX_GC, $ARGV[1] or die $!;
open my $compOut, '>', $ARGV[2] or die $!;
open my $dupeOut, '>', $ARGV[3] or die $!;
open my $netOut, '>', $ARGV[4] or die $!;

while(<IDX_GC>)
{
	chomp;
	my @sps = split(/\t/);
	my @sps2 = split(/:/, $sps[0]);
	$sps[0] = $sps2[0];
	#print $sps[0] . "\n";
	$sps[3] = ($sps[2] * 150) / $sps[1];
	
	@{$idxHash{$sps[0]}} = @sps;
}

close IDX_GC;

my @keys = keys %hitHash;

foreach(@keys)
{
	my $bus = $_;
	my @hits = @{$hitHash{$bus}};
	
	if(scalar @hits == 1)
	{
		my @idx = @{$idxHash{$hits[0]}};
		
		print {$compOut} join("\t", @idx) . "\n";
	}
	elsif(scalar @hits > 1)
	{
		foreach(@hits)
		{
			my $hit = $_;
			
			my @idx = @{$idxHash{$hit}};
			print {$dupeOut} join("\t", @idx) . "\n";
		}
		
		my @x0; my @x1; my @y0; my @y1;
		
		for(my $i = 0; $i < scalar @hits - 1; $i++)
		{
			my @idx1 = @{$idxHash{$hits[$i]}};
			
			if($idx1[1] > 5000)
			{
				for(my $j = $i + 1; $j < scalar @hits; $j++)
				{
					my @idx2 = @{$idxHash{$hits[$j]}};
					
					if($idx2[1] > 5000)
					{
						push(@x0, $idx1[4]); push(@x1, $idx2[4]);
						push(@y0, $idx1[3]); push(@y1, $idx2[3]);
					}
				}
			}
		}
		
		for(my $i = 0; $i < scalar @x0; $i++)
		{
			print {$netOut} $x0[$i] . "\t" . $y0[$i] . "\t" . $x1[$i] . "\t" . $y1[$i] . "\n";
		}
	}
}
close $compOut;
close $dupeOut;
close $netOut;