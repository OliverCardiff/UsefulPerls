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
		push(@{$hitHash{$sps[2]}}, $sps[0]);
	}
	
}

close BUSCO;

my %scaffCov;
my %scaffGC;
my %scaffSize;

open IDX_GC, $ARGV[1] or die $!;

while(<IDX_GC>)
{
	chomp;
	my @sps = split(/\t/);
	my @sps2 = split(/:/, $sps[0]);
	$sps[0] = $sps2[0];

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
my $testNo = 0;

my %OptiA;
my %OptiB;

my %OptiALevels;
my %OptiBLevels;

my %OptiACrit;
my %OptiBCrit;

my @blankLvls = (0,0,0,0,0);
my @blankCrit = (0,0,0,0,0,0,0);

for(my $i = 0.5; $i < 3.2; $i += 0.2)
{
	$OptiA{$i} = 0;
	$OptiB{$i} = 0;	
	
	@{$OptiALevels{$i}} = @blankLvls;
	@{$OptiBLevels{$i}} = @blankLvls;
	
	@{$OptiACrit{$i}} = @blankCrit;
	@{$OptiBCrit{$i}} = @blankCrit;
}


for(my $covPos = 0.75; $covPos < 1.26; $covPos += 0.05)
{
	for(my $covNeg = 0.75; $covNeg < 1.26; $covNeg += 0.05)
	{
		for(my $GCPos = 0.75; $GCPos < 1.26; $GCPos += 0.05)
		{
			for(my $GCNeg = 0.75; $GCNeg < 1.26; $GCNeg += 0.05)
			{
				for(my $ovalGC = 39.2; $ovalGC < 40.9; $ovalGC += 0.2)
				{
					for(my $ovalCov = 92; $ovalCov < 103; $ovalCov++)
					{
						for(my $limit = 16; $limit < 25; $limit++)
						{
							my @chosen;
							my $sSize = 0;
							
							foreach(@scaffs)
							{
								my $sc = $_;
								
								my $res = CheckOvalInclude($covPos, $covNeg, $GCPos, $GCNeg, $ovalGC, $ovalCov, $limit, $scaffCov{$sc}, $scaffGC{$sc});
								
								if($res == 1)
								{
									push(@chosen, $sc);
									$sSize += $scaffSize{$sc};
								}
							}
							
							my @levels = TestSetDC(\@chosen);
							push(@levels, $sSize);
							push(@levels, $sSize);
							$levels[3] -= 500000000;
							if($levels[3] < 0) 
							{
								$levels[3] = 0;
							}
							$levels[3] /= 3000000;
							my @values = ($covPos, $covNeg, $GCPos, $GCNeg, $ovalGC, $ovalCov, $limit);
							
							##Selecting the Optimisation terms
							
							for(my $i = 0.5; $i < 3.2; $i += 0.2)
							{
								if($levels[0] - ($levels[1] * $i) - $levels[3] > $OptiA{$i})
								{
									$OptiA{$i} = $levels[0] - ($levels[1] * $i) - $levels[3];
									@{$OptiALevels{$i}} = @levels;
									@{$OptiACrit{$i}} = @values;
								}
								
								if($levels[0] - ($levels[2] * $i) - $levels[3] > $OptiB{$i})
								{
									$OptiB{$i} = $levels[0] - ($levels[2] * $i) - $levels[3];
									@{$OptiBLevels{$i}} = @levels;
									@{$OptiBCrit{$i}} = @values;
								}
							}

							
							$testNo++;
							
							if($testNo % 100 == 0)
							{
								my $tval = 1.5;
								
								print "Tests Ran: " . $testNo . "\n";
								print "Complete - Single-Dupe Optimisation Progress: " . $OptiA{$tval} . "\n";
								print "Results:\t" . join("\t", @{$OptiALevels{$tval}}) . "\n";
								print "Terms:  \t" . join("\t", @{$OptiACrit{$tval}}) . "\n";
								print "Complete - Total-Dupe Optimisation Progress: " . $OptiB{$tval} . "\n";
								print "Results:\t" . join("\t", @{$OptiBLevels{$tval}}) . "\n";
								print "Terms:  \t" . join("\t", @{$OptiBCrit{$tval}}) . "\n";

								print "\nNow:    \t" . join("\t", @values) . "\n\n";
							}
						}
					}
				}
			}
		}
	}
}

open my $fout, '>', $ARGV[2] or die $!;

for(my $i = 0.5; $i < 3.2; $i += 0.2)
{
	print {$fout} $i . "\t" . $OptiA{$i} . "\t" . join("\t", @{$OptiALevels{$i}}) . "\t" . join("\t", @{$OptiACrit{$i}}) . "\n";
}

close $fout;

open my $fout2, '>', $ARGV[3] or die $!;

for(my $i = 0.5; $i < 3.2; $i += 0.2)
{
	print {$fout2} $i . "\t" . $OptiB{$i} . "\t" . join("\t", @{$OptiBLevels{$i}}) . "\t" . join("\t", @{$OptiBCrit{$i}}) . "\n";
}

close $fout2;

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

sub TestSetDC
{
	my %buscoHash;
	my %dupeHash;
	my $totDupe = 0;
	my $totUni = 0;
	
	my @array = @{$_[0]};
	
	foreach(@array)
	{
		my $sc = $_;
		my @buscs = ();
		if(exists $hitHash{$sc})
		{
			@buscs = @{$hitHash{$sc}};
		}
		if((scalar @buscs) > 0)
		{
			foreach(@buscs)
			{
				my $busc = $_;
				
				if(exists $buscoHash{$busc})
				{
					$buscoHash{$busc}++;
					$dupeHash{$busc}++;
					$totDupe++;
				}
				else
				{
					$buscoHash{$busc} = 1;
					$totUni++;
				}
			}
		}
	}
	
	my @keysA = keys %dupeHash;
	my $uniDupe = scalar @keysA;
	my @return = ($totUni, $uniDupe, $totDupe);
	return @return;
}
