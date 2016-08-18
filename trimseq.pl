use strict;
use warnings;
use POSIX;
use List::Util qw(sum);

sub median {
  sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), ceil( $#_/2 ) ] )/2;
}

my $cutSize = $ARGV[2];

open CONTIGS, $ARGV[0] or die $!;

open my $fout, '>', $ARGV[1] or die $!;
#open my $scrapOut, '>', $ARGV[2] or die $!;

my $oldID = "";
my $seq = "";
my $size = 0;
my %namHash;
my @contigs;
my %namHash2;

while(<CONTIGS>)
{
	chomp;
	my $subs = substr $_, 0 , 1;
	
	if($subs eq ">")
	{
		if(!($oldID eq ""))
		{
			if (length($seq) > $cutSize)
			{
				$namHash{$oldID} = length($seq);
				$size += length($seq);
				push(@contigs, $namHash{$oldID});
			}
			elsif(length($seq) > 500)
			{
				$namHash2{$oldID} = 1;
			}
		}
		$seq = "";
		$oldID = $_;
	}
	else
	{
		$seq = $seq . $_;
	}
}

close CONTIGS;



open CONTIGS, $ARGV[0] or die $!;

my $subPrint = 0;
my $printBool = 0;
while(<CONTIGS>)
{
	chomp;
	my $subs = substr $_, 0 , 1;
	my $line = $_;
	
	if($subs eq ">")
	{
		$seq = "";
		$oldID = $_;
		
		if(exists $namHash{$oldID})
		{
			$printBool = 1;
			my @sps = split(/\|/);
			
			$line = $sps[0] . "_length:" . $namHash{$oldID};
		}
		elsif(exists $namHash2{$oldID})
		{
			$subPrint = 1;
			$printBool = 0;
		}
		else
		{
			$subPrint = 0;
			$printBool = 0;
		}
	}
	
	if($printBool == 1)
	{
		print {$fout} $line . "\n";
	}
	elsif($subPrint == 1)
	{
#		print {$scrapOut} $line . "\n";
	}
}

close $fout;
#close $scrapOut;
close CONTIGS;


print "Total > $cutSize: " . scalar (keys %namHash) . "\n"; 
print "N50: " . median(@contigs) . "\n";
print "Genome Size: " . $size . "\n";
