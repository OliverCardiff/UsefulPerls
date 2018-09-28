use strict;
use warnings;

open CDS, $ARGV[0] or die $!;
open my $annot, '>', $ARGV[1] or die $!;

my $size = 0;
my $pos = 0;
my $scaff = 0;

while(<CDS>)
{
	chomp;
	
	my $subs = substr $_, 0, 1;
	
	if($subs eq ">")
	{
		my @sps = split(/ /);
		my @sps2 = split(/:/, $sps[4]);
		
		my @scf = split(/\|/, $sps2[0]);
		
		my @numsA = split(/\(/, $sps2[1]);
		my @numsB = split(/-/, $numsA[0]);
		
		if($numsB[0] > $numsB[1])
		{
			my $temp = $numsB[0];
			$numsB[0] = $numsB[1];
			$numsB[1] = $temp;
		}
		
		$size = $numsB[1] - $numsB[0];
		$pos = $numsB[0] + ($size/2);
		$scaff = $scf[0];
	} else
	{
		my $count = ($_ =~ tr/X//);
		if($count < 20)
		{
			print {$annot} $scaff . "\t" . $size . "\t" . $pos . "\t1\n";
		}
	}
}

close $annot;
close CDS;