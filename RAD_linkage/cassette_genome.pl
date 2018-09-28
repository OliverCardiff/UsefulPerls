use strict;
use warnings;

my %spanHash;
my %blockFiles;
my $stub = $ARGV[2];

open NAMVEC, $ARGV[0] or die $!;

while(<NAMVEC>)
{
	chomp;

	my @sps = split(/\t/);

	$spanHash{$sps[1]}{$sps[0]} = 1;
}

close NAMVEC;

my @block_keys = keys %spanHash;

foreach(@block_keys)
{
	my $block = $_;

	open my $fout, '>', "LR_" . $stub . "_" . $block . ".fa" or die $!;

	$blockFiles{$block} = $fout;
}

open GENOME, $ARGV[1] or die $!;

my $printr;
my $oldID = "";
my $switch = 0;

while(<GENOME>)
{
	chomp;

	my $subs = substr $_, 0, 1;
	my $line = $_;

	if($subs eq ">")
	{
		my $nam = $_;
		$nam =~ tr/>//d;
		
		my @blk = keys %spanHash;

		$switch = 0;

		foreach(@blk)
		{
			if(exists $spanHash{$_}{$nam})
			{
				$printr = $blockFiles{$_};
				$switch = 1;
			}
		}
		
	}
	
	if($switch == 1)
	{
		print {$printr} $line . "\n";
	}
}

close GENOME;

my @blk = keys %spanHash;
foreach(@blk)
{
	close $blockFiles{$_};
}

