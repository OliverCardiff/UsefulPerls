# Usage
# perl contig_select.pl <list.txt> <genome.fasta> <listcontigs.fasta> <nonlistcontigs.fasta>

use strict;
use warnings;

my %scaffHash;

open LIST, $ARGV[0] or die $!;

while(<LIST>)
{
	chomp;
	my $scaff = $_;
	$scaffHash{$scaff} = 1;
}

close LIST;

open GENOME, $ARGV[1] or die $!;
open my $listOut, '>', $ARGV[2] or die $!;
open my $nonOut, '>', $ARGV[3] or die $!;

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
		
		if(exists $scaffHash{$nam})
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
		print {$listOut} $line . "\n";
	}
	else
	{
		print {$nonOut} $line . "\n";
	}
}

close GENOME;
close $listOut;
close $nonOut;
