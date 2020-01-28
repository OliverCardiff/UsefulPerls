#usages matchNNN.pl <reflistNR.txt> <reflistNG.txt> <hollow.fa> <gapfilled.fa> <reference.fa>
use strict;
use warnings;

my %nnHash;
my %gapHash;
my %refHash;
my %alnNG;
my %alnNR;

my %scSts;
my %scEds;

open F1, $ARGV[2] or die $!;
my $oldID = "";
my $seq = "";

my @ids;

while(<F1>)
{
	chomp;

	my $subs = substr $_, 0, 1;
	
	if($subs eq ">")
	{
		if($oldID ne "")
		{
			$nnHash{$oldID} = $seq;
			push(@ids, $oldID);
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
push(@ids, $oldID);
$nnHash{$oldID} = $seq;
close F1;

open F2, $ARGV[3] or die $!;
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
			$gapHash{$oldID} = $seq;
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
$gapHash{$oldID} = $seq;
close F2;

open F3, $ARGV[4] or die $!;
$oldID = "";
$seq = "";
while(<F3>)
{
	chomp;

	my $subs = substr $_, 0, 1;
	
	if($subs eq ">")
	{
		if($oldID ne "")
		{
			$refHash{$oldID} = $seq;
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
$refHash{$oldID} = $seq;
close F3;

open REF, $ARGV[0] or die $!;

while(<REF>)
{
	chomp;
	my @sps = split(/\t/);
	if($sps[2] > 1)
	{
		$alnNR{$sps[0]} = $sps[1];
	}
}

close REF;

open REF2, $ARGV[1] or die $!;

while(<REF2>)
{
	chomp;
	my @sps = split(/\t/);
	if($sps[2] > 1)
	{
		$alnNG{$sps[0]} = $sps[1];
	}
}

close REF2;

my @keyset = keys %nnHash;

foreach(@keyset)
{
	my $sc = $_;
	#print $sc . "\n";
	my $lim = length($nnHash{$sc});
	my $seq = $nnHash{$sc};
	my $track = 0;
	my $isN = 0;

	my $cnt = 0;
	my $cnt2 = 0;

	for(my $i = 0; $i < $lim; $i++)
	{
		my $ch = substr $seq, $i , 1;
		
		if($ch ne "N")
		{
			if($isN == 1)
			{
				push(@{$scEds{$sc}}, $i);
				$cnt2++;
				$isN = 0;
			}
		} else {
			
			if($isN == 0)
			{
				push(@{$scSts{$sc}}, $i);
				$cnt++;
				$isN = 1;
			}					
		}
		$track = $i;
	}
	if($isN == 1)
	{
		push(@{$scEds{$sc}}, $track);
		$cnt2++;
	}


	print $sc . "\t" . "sts: " . $cnt . "\t" . "eds: " . $cnt2 . "\n";
}

@keyset = keys %nnHash;

foreach(@keyset)
{
	my $sc = $_;
	if(exists $alnNR{$sc} && exists $alnNG{$sc})
	{
		my $extNR = $alnNR{$sc};
		my $extNG = $alnNG{$sc};
		my $lim = length($nnHash{$sc});
		my $seqR = $refHash{$sc};
		my $seqG = $gapHash{$sc};

		my @sts = @{$scSts{$sc}};
		my @eds = @{$scEds{$sc}};
		my $ln = scalar @sts;

		for(my $i = 0; $i < $ln; $i++)
		{
			my $st = $sts[$i] - 5000;
			my $ed = $eds[$i] + 5000;
			my $siz = $ed-$st;
			
			if(($st - $extNG) > 0 && ($st - $extNR) > 0)
			{
				my $subG = substr $seqG, ($st - $extNG), ($siz);
				my $subR = substr $seqR, ($st - $extNR), ($siz);

				my $cnt = ($subG =~ tr/N//);
				my $l1 = length($subG);
				my $l2 = length($subR);
		
				if($cnt < ($siz * 0.1) && $siz < 16000 && $siz > 12000 && $l1 == $l2)
				{
					open my $fout, '>', "allele_frags/$sc" . "_$i" . ".fa";
			
					print {$fout} ">Gapped_$st" . "_$ed\n";
					print {$fout} $subG . "\n";
					print {$fout} ">Reference_$st" . "_$ed\n";	
					print {$fout} $subR . "\n";	
						
					close $fout;
				}
			}
		}
	}
}
