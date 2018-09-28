# || ARGV[0] IDXstats || ARGV[1] Genome || ARGV[2] GenomeOut || ARGV[3] StatsOut ||

use strict;
use warnings;

my %scaffCov;
my %scaffGC;
my %scaffSize;

open IDX_GC, $ARGV[0] or die $!;

while(<IDX_GC>)
{
	chomp;
	my @sps = split(/\t/);
	#my @sps2 = split(/:/, $sps[0]);
	#$sps[0] = $sps2[0];

	$sps[3] = ($sps[2] * 150) / $sps[1];
	
	if($sps[3] > 55 && $sps[3] < 130 && $sps[4] > 36 && $sps[4] < 44)
	{
		$scaffCov{$sps[0]} = $sps[3];
		$scaffGC{$sps[0]} = $sps[4];
		$scaffSize{$sps[0]} = $sps[1];
	}	
}

close IDX_GC;

my @scaffs = keys %scaffCov;
my $scCount = scalar @scaffs;
print $scCount . "\n";

#0.8	0.75	0.85	0.75	40	98	24
my $covPos = 0.8; my $covNeg = 0.75; my $GCPos = 0.85;
my $GCNeg = 0.75; my $ovalGC = 40; my $ovalCov = 98;
my $limit = 24;

my %chosen;
my $sSize = 0;

foreach(@scaffs)
{
	my $sc = $_;
	
	my $res = CheckOvalInclude($covPos, $covNeg, $GCPos, $GCNeg, $ovalGC, $ovalCov, $limit, $scaffCov{$sc}, $scaffGC{$sc});
	
	if($res == 1)
	{
		$chosen{$sc} = 1;
		#print $sc . "\n";
		$sSize += $scaffSize{$sc};
	}
}

open GENOME, $ARGV[1] or die $!;
open my $fout, '>', $ARGV[2] or die $!;
my $switch = 0;

while(<GENOME>)
{
	chomp;
	my $line = $_;
	
	my $subs = substr $_, 0, 1;
	
	if($subs eq ">")
	{
		my $nam = $line;
		$nam =~ tr/>//d;

		if(exists $chosen{$nam})
		{
			$switch = 1;
		}
		else
		{
			$switch = 0;
		}
	}
	
	if($switch == 1)
	{
		print {$fout} $line . "\n";
	}
}

close $fout;
close GENOME;

open IDX, $ARGV[0] or die $!;
open my $fout2, '>', $ARGV[3] or die $!;
my $switch = 0;

while(<IDX>)
{
	chomp;
	my $line = $_;
	
	my @sps = split(/\t/);
	
	if(exists $chosen{$sps[0]})
	{
		print {$fout2} $line;
	}
}

close $fout2;
close IDX;


sub CheckOvalInclude
{
	my ($covPos, $covNeg, $GCPos, $GCNeg, $ovalGC, $ovalCov, $limit, $testCov, $testGC) = (@_);
	
	my $dCov = $testCov - $ovalCov;
	my $dGC = $testGC - $ovalGC;
	
	$dGC *= 10;
	
	if($dCov > 0)
	{
		$dCov *= $covPos;
	}
	if($dCov < 0)
	{
		$dCov *= $covNeg;
	}
	if($dGC > 0)
	{
		$dGC *= $GCPos;
	}
	if($dGC < 0)
	{
		$dGC *= $GCNeg;
	}
	
	$dGC *= $dGC;
	$dCov *= $dCov;
	
	if(($dGC + $dCov) < ($limit*$limit))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
