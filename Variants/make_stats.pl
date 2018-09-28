use strict;
use warnings;

my %lengthHash;

open FAI, $ARGV[0] or die $!;

while(<FAI>)
{
	chomp;
	my @sps = split(/\t/);
	
	if(scalar @sps > 1)
	{
		$lengthHash{$sps[0]} = $sps[1];
	}
}

close FAI;

my %wigHash;
my %countHash;

my $fiftyPlus = 0;
my %fiftyHash;
my $oldID  = "";

open VCF, $ARGV[1] or die $!;

while(<VCF>)
{
	chomp;
	
	if($. > 54)
	{
		my @sps = split(/\t/);
		my $len = length($sps[3]);
		my @sps2 = split(/_/, $sps[0]);
		
		if($oldID eq $sps[0])
		{
			while($sps[1] > $fiftyPlus)
			{
				$fiftyPlus += 50;
				$fiftyHash{$sps2[0]}{$fiftyPlus} = 0;
			}
			
			$fiftyHash{$sps2[0]}{$fiftyPlus} += $len;
		} else {
			$fiftyPlus = 50;
			$fiftyHash{$sps2[0]}{$fiftyPlus} = 0;
			print "Did " . $sps2[0] . "\n";
			while($sps[1] > $fiftyPlus)
			{
				$fiftyPlus += 50;
				$fiftyHash{$sps2[0]}{$fiftyPlus} = 0;
			}
			
			$fiftyHash{$sps2[0]}{$fiftyPlus} += $len;
		}
		
		if(exists $countHash{$sps[0]})
		{
			$countHash{$sps[0]} += $len;
		} else {
			$countHash{$sps[0]} = $len;
		}
		
		$oldID = $sps[0];
	}
	
}


close VCF;

my @keys = keys %fiftyHash;

open my $fout, '>', $ARGV[3] or die $!;

foreach(@keys)
{
	my $scaf = $_;
	
	my $fif = 50;
	
	while(exists $fiftyHash{$scaf}{$fif})
	{	
		print {$fout} $scaf . "\t" . $fiftyHash{$scaf}{$fif} . "\n";
		$fif += 50;
	}
}

close $fout;

open my $countOut, '>', $ARGV[2] or die $!;
my @keys2 = keys %countHash;

foreach(@keys2)
{
	my $scaf = $_;
	my $len = $lengthHash{$scaf};
	my $cnt = $countHash{$scaf};
	
	my $score = ($cnt / $len) * 100;
	
	print {$countOut} $scaf . "\t" . $len . "\t" . $cnt . "\t" . $score . "\n";
}

close $countOut;