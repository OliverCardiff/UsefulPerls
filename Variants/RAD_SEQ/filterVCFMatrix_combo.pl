use strict;
use warnings;

my $pHigh = 2.5;
my $pLow = 0.9;

my $iHigh = 0.02;
my $iLow = 0.01;

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

open my $highOut, '>', $ARGV[2] or die $!;
open my $lowOut, '>', $ARGV[3] or die $!;
 
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

			if($poly > $pHigh && $intra > $iHigh)
			{
				print {$highOut} $line . "\n";
			}

			if($poly < $pLow && $intra < $iLow)
			{
				print {$lowOut} $line . "\n";
			}
		}
	}
	else {
		
		print {$highOut} $line . "\n";
		print {$lowOut} $line . "\n";	
	}
}

close VCF;
close $highOut;
close $lowOut;
