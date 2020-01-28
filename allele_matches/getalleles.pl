use strict;
use warnings;

my %oneHash;
my %twoHash;

open F1, $ARGV[0] or die $!;
my $oldID = "";
my $seq = "";
while(<F1>)
{
	chomp;

	my $subs = substr $_, 0, 1;
	
	if($subs eq ">")
	{
		if($oldID ne "")
		{
			$oneHash{$oldID} = $seq;
		}

		my @sps = split(/\./);
		$oldID = $sps[0];
		$oldID =~ tr/>//d;
		#print $oldID . "\n";
		$seq = "";
	}
	else {
		$seq .= $_;	
	}

}
$oneHash{$oldID} = $seq;
close F1;

open F2, $ARGV[1] or die $!;
$oldID = "";
$seq = "";
while(<F2>)
{
	chomp;

	my $subs = substr $_, 0, 1;
	
	if($subs eq ">")
	{
		if($oldID ne "")
		{
			$twoHash{$oldID} = $seq;
		}

		my @sps = split(/\./);
		$oldID = $sps[0];
		$oldID =~ tr/>//d;
		#print $oldID . "\n";
		$seq = "";
	}
	else {
		$seq .= $_;	
	}

}
$oneHash{$oldID} = $seq;
close F2;

my @keyset = keys %oneHash;

foreach(@keyset)
{
	my $sc = $_;

	open my $fout1, '>', "NR1/$sc.fa" or die $!;
	open my $fout2, '>', "NR2/$sc.fa" or die $!;
	print "testing: " . $sc . "\n";
	print {$fout1} ">" . $sc . "_1\n";
	print {$fout1} $oneHash{$sc} . "\n";
	print {$fout2} ">" . $sc . "_2\n";
	print {$fout2} $twoHash{$sc} . "\n";

	close $fout1;
	close $fout2;
}


