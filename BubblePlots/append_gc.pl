use strict;
use warnings;

my %gcHash;

open GENOME, $ARGV[0] or die $!;

my $concat = "";
my $oldID = "";
my $bases = 0;
my $thresh = 10000000;

while(<GENOME>)
{
	chomp;
	chomp;
	my $line = $_;
	my $subs = substr $_, 0, 1;
	
	if($subs eq ">")
	{
		if($bases > $thresh)
		{
			$thresh += 10000000;
			print $bases . " bases read so far\n";
		}
		my $nam = $line;
		$nam =~ tr/>//d;
		
		$gcHash{$nam} = 0;
		
		if($oldID ne "")
		{
			my $len = length($concat);
			my $at = 0; my $gc = 0;
			for(my $i = 0; $i < $len; $i++)
			{
				my $ch = substr $concat, $i, 1;
				if($ch eq "G" || $ch eq "C" || $ch eq "g" || $ch eq "c")
				{
					$gc++;
				}
				elsif($ch eq "A" || $ch eq "T" || $ch eq "a" || $ch eq "t")
				{
					$at++;
				}
				$bases++;
			}
			
			my $ratio = $gc/($at + $gc) * 100;
			$gcHash{$oldID} = $ratio;
		}
		
		$oldID = $nam;
		$concat = "";
	}
	else
	{
		$concat .= $line;
	}
	
}

my $len = length($concat);
my $at = 0; my $gc = 0;
for(my $i = 0; $i < $len; $i++)
{
	my $ch = substr $concat, $i, 1;
	if($ch eq "G" || $ch eq "C" || $ch eq "g" || $ch eq "c")
	{
		$gc++;
	}
	elsif($ch eq "A" || $ch eq "T" || $ch eq "a" || $ch eq "t")
	{
		$at++;
	}
}

my $ratio = $gc/($at + $gc) * 100;
$gcHash{$oldID} = $ratio;

close GENOME;

open IDX, $ARGV[1] or die $!;
open my $fout, '>', $ARGV[2] or die $!;

while(<IDX>)
{
	chomp;
	my @sps = split(/\t/);
	
	if(exists $gcHash{$sps[0]})
	{
		print {$fout} join("\t", @sps) . "\t" . $gcHash{$sps[0]} . "\n";
	}
}

close IDX;
close $fout;