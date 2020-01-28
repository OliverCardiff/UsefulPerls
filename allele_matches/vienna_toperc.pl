use strict;
use warnings;

my $seq1 = "";
my $seq2 = "";

open VIE, $ARGV[0] or die $!;
my $cnt = 0;
while(<VIE>)
{
	chomp;
	if($. == 2)
	{
		$seq1 = $_;
	}
	if($. == 4)
	{
		$seq2 = $_;
	}
	
}

close VIE;

$seq1 =~ tr/atgc/ATGC/;
$seq2 =~ tr/atgc/ATGC/;

my $len = length($seq1);
my $started = 0;
my $innerCnt = 0;
my $matchCnt = 0;
my $misMatches = 0;
my $totIndel = 0;
my $gapAccu = 0;

my @matchList;
my @errList;
my @sc1List;
my @sc2List;
my @inds;

my $ind = 0;
my @tInds = ();
my $maxInd = 0;

my $trSeq1 = "";
my $trSeq2 = "";


for(my $i = 0; $i < $len; $i++)
{
        my $ch1 = substr $seq1, $i, 1;
        my $ch2 = substr $seq2, $i, 1;

	if($started == 0)
	{
		if($ch1 ne "-" && $ch2 ne "-")
		{
			$started = 1;
		}
	}
	else
	{
		if($ch1 ne "-")
		{
			push(@sc1List, $ind);
			$trSeq1 .= $ch1;
		}
		if($ch2 ne "-")
		{
			push(@sc2List, $ind);
			$trSeq2 .= $ch2;
		}
		if($ch1 eq $ch2)
		{
			if($gapAccu > 0)
			{
				$totIndel += $gapAccu;
				$innerCnt += $gapAccu;
				
				for(my $j = 0; $j < $gapAccu; $j++)
				{
					push(@matchList, 1);
					push(@errList, 0);
				}
				foreach(@tInds)
                                {
                                	push(@inds, $_);
                                }
				$gapAccu = 0;
                                @tInds = ();
			}
			$maxInd = $ind;
			$matchCnt++;
                        push(@matchList, 1);
                        push(@errList, 1);
                        push(@inds, $ind);
                        $innerCnt++;
		}
		else
		{
			if($ch1 ne "-" && $ch2 ne "-")
         	        {
				if($gapAccu > 0)
                        	{
                                	$totIndel += $gapAccu;
                                	$innerCnt += $gapAccu;
                                	
					for(my $j = 0; $j < $gapAccu; $j++)
                                	{
                                        	push(@matchList, 1);
						push(@errList, 0);
                               		}
					foreach(@tInds)
					{
						push(@inds, $_);
					}
					$gapAccu = 0;
					@tInds = ();
                        	}
				$misMatches++;
                                push(@matchList, 0);
                                push(@errList, 0);
                                push(@inds, $ind);
                                $innerCnt++;
				$maxInd = $ind;
			}
			else
			{
				$gapAccu++;
				push(@tInds, $ind);
			}
		}
		$ind++;
	}
}

my $totalErr = $totIndel + $misMatches;
my $perc = ($matchCnt / $innerCnt) * 100;

system("mkdir vienna");

open my $sumOut, '>', "vienna/summary.txt" or die $!;

my $sc1Len = scalar @sc1List;
my $sc2Len = scalar @sc2List;

print {$sumOut} $innerCnt . "\t" . $matchCnt . "\t" . $perc . "\t" . $totIndel . "\t" . $misMatches . "\t" . $sc1Len . "\t" . $sc2Len  ."\n";

close $sumOut;
#print "Identity is $perc\n";
#print "Aligned Size is $innerCnt\n";
#print "There were $totIndel bases of indel seq\n";
#print "There were $misMatches nucleotide substitutions\n";

open my $sc1Out, '>', "vienna/sc1.txt" or die $!;
my $sqCnt1 = 0;
foreach(@sc1List)
{
	if($_ < $maxInd)
	{
		print {$sc1Out} $_ . "\n";
		$sqCnt1++;
	}
}

close $sc1Out;

open my $sc2Out, '>', "vienna/sc2.txt" or die $!;
my $sqCnt2 = 0;
foreach(@sc2List)
{
	if($_ < $maxInd)
	{
        	print {$sc2Out} $_ . "\n";
		$sqCnt2++;
	}
}

close $sc2Out;

open my $sqfrag1, '>', "vienna/frag1.fa" or die $!;

print {$sqfrag1} ">Frag1\n";
my $toprint = substr $trSeq1, 0, $sqCnt1;
print {$sqfrag1} $toprint . "\n";

close $sqfrag1;

open my $sqfrag2, '>', "vienna/frag2.fa" or die $!;

print {$sqfrag2} ">Frag2\n";
$toprint = substr $trSeq2, 0, $sqCnt2;
print {$sqfrag2} $toprint . "\n";

close $sqfrag2;

open my $identOut, '>', "vienna/idents.txt" or die $!;

my $len1 = scalar @matchList;

print {$identOut} "Index\tSubstitutes\tTotalErr\n";
for(my $i = 0; $i < $len1; $i++)
{
	print {$identOut} $inds[$i] . "\t" . $matchList[$i] . "\t" . $errList[$i] . "\n";
}

close $identOut;
