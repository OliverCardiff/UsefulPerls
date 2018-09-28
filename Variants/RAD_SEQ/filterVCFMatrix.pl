use strict;
use warnings;

my $pHigh = 6;
my $pLow = 1;

my $iHigh = 0.49;
my $iLow = 0.02;

my %polyHash;
my %intraHash;

open MATRIX, $ARGV[0] or die $!;

while(<MATRIX>)
{
	chomp;
	my @sps = split(/\t/);

	my $nam = $sps[0] . "_" . $sps[1];

	$polyHash{$nam} = $sps[2];
	$intraHash{$nam} = $sps[3];
}

close MATRIX;

open my $pOut, '>', $ARGV[2] or die $!;
open my $iOut, '>', $ARGV[3] or die $!;
 
open VCF, $ARGV[1] or die $!;

while(<VCF>)
{
	chomp;
	my $line = $_;
	my $subs = substr $_, 0, 1;

	if($subs ne "#")
	{
		my @sps = split(/\t/);
		my @sps2 = split(/\|/, $sps[0]);

		my $nam = $sps2[0] . "_" . $sps[1];

		if(exists($polyHash{$nam}))
		{
			my $poly = $polyHash{$nam};
			my $intra = $intraHash{$nam};

			if($poly > $pLow && $poly < $pHigh)
			{
				print {$pOut} $line . "\n";
			}

			if($intra > $iLow && $intra < $iHigh)
			{
				print {$iOut} $line . "\n";
			}
		}
	}
	else {
		
		print {$pOut} $line . "\n";
		print {$iOut} $line . "\n";	
	}
}

close VCF;
close $pOut;
close $iOut;
