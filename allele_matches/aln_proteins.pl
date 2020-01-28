use strict;
use warnings;

my $folds = $ARGV[1];
my $stub = $ARGV[0];

my $success = 1;

for(my $i = 1; $i <= $folds; $i++)
{
	if($success == 0)
	{
		chdir("../../../");
	}
	$success = 0;
	
	chdir($stub . "/Pairs_" . $i . "/vienna");
	
	my @fr1s = ();
	my @fr2s = ();
	my %fr1Keeps;
	my %fr2Keeps;

	system("pwd");	
	print "opening fr1.cds\n";
	open FR1a, "fr1.cds";
	while(<FR1a>)
	{
		chomp;
		push(@fr1s, $_);
	}
	close FR1a;
	print "opening fr1.cds\n";
	open FR2a, "fr2.cds";
	while(<FR2a>)
	{
		chomp;
		push(@fr2s, $_);
	}
	close FR2a;
	
	my $fr1Cnt = scalar @fr1s;
	my $fr2Cnt = scalar @fr2s;
	
	if($fr1Cnt < 1 || $fr2Cnt < 1)
	{
		print "Not enough ORFS\n";
		next;
	}
	
	if($fr1Cnt > 1)
	{
		for(my $j = 0; $j < $fr1Cnt; $j++)
		{
			my $keep = 1;
			for(my $k = 0; $k < $fr1Cnt; $k++)
			{
				if($j != $k)
				{
					my @spJ = split(/\t/, $fr1s[$j]);
					my @spK = split(/\t/, $fr1s[$k]);
					
					if(($spJ[0] > $spK[0] && $spJ[0] < $spK[1]) || ($spJ[1] < $spK[1] && $spJ[1] > $spK[0]))
					{
						my $lenJ = $spJ[1] - $spJ[0];
						my $lenK = $spK[1] - $spK[0];
						if($lenJ < $lenK)
						{
							$keep = 0;
						}
					}
				}
			}
			my $jind = $j + 1;
			$fr1Keeps{$jind} = $keep;
		}
	}
	else
	{
		$fr1Keeps{1} = 1;
	}
	
	if($fr2Cnt > 1)
	{
		for(my $j = 0; $j < $fr2Cnt; $j++)
		{
			my $keep = 1;
			for(my $k = 0; $k < $fr2Cnt; $k++)
			{
				if($j != $k)
				{
					my @spJ = split(/\t/, $fr2s[$j]);
					my @spK = split(/\t/, $fr2s[$k]);
					
					if(($spJ[0] > $spK[0] && $spJ[0] < $spK[1]) || ($spJ[1] < $spK[1] && $spJ[1] > $spK[0]))
					{
						my $lenJ = $spJ[1] - $spJ[0];
						my $lenK = $spK[1] - $spK[0];
						if($lenJ < $lenK)
						{
							$keep = 0;
						}
					}
				}
			}
			my $jind = $j + 1;
			$fr2Keeps{$jind} = $keep;
		}
	}
	else
	{
		$fr2Keeps{1} = 1;
	}
	
	my $cnt1 = 1;
	
	my @fr1Keys = keys %fr1Keeps;
	my @fr2Keys = keys %fr2Keeps;
	
	print "opening fr1x.cds\n";;
	open my $cds1, '>', "fr1x.cds";
	
	foreach(@fr1Keys)
	{
		my $fr1I = $_;
		if($fr1Keeps{$_} == 1)
		{
			print {$cds1} $fr1s[$fr1I-1] . "\n";
		}
	}
	
	close $cds1;
	print "opening fr2x.cds\n";
	open my $cds2, '>', "fr2x.cds";
	
	foreach(@fr2Keys)
	{
		my $fr2I = $_;
		if($fr2Keeps{$_} == 1)
		{
			print {$cds2} $fr2s[$fr2I-1] . "\n";
		}
	}
	
	close $cds2;

	my $printedPeps1 = 0;
	my $printedPeps2 = 0;
	
	print "opening allpeps.pep\n";
	open my $allpeps, '>', "allpeps.pep";
	print "opening frag1.fa.transdecoder_dir/longest_orfs.pep\n";
	open PEPS1, "frag1.fa.transdecoder_dir/longest_orfs.pep" or next;
	my $switch = 0;
	my $cntr = 1;

	my $pepNam1 = "";
	my $pepNam2 = "";
	
	while(<PEPS1>)
	{
		chomp;
		my $line = $_;
		my $subs = substr $_, 0, 1;
		if($subs eq ">")
		{
			if($fr1Keeps{$cntr} == 1)
			{
				$switch = 1;
				my @sps = split(/\s+/);
				$pepNam1 = $sps[0];
				$pepNam1 =~ tr/>//d;
				$printedPeps1++;
			}
			else
			{
				$switch = 0;
			}
			$cntr++;
		}
		if($switch == 1)
		{
			print {$allpeps} $line . "\n";
		}
	}
	close PEPS1;
	print "opening frag2.fa.transdecoder_dir/longest_orfs.pep\n";
	open PEPS2, "frag2.fa.transdecoder_dir/longest_orfs.pep" or next;
	$cntr = 1;
	
	while(<PEPS2>)
	{
		chomp;
		my $line = $_;
		my $subs = substr $_, 0, 1;
		if($subs eq ">")
		{
			if($fr2Keeps{$cntr} == 1)
			{
				$switch = 1;
				my @sps = split(/\s+/);
				$pepNam2 = $sps[0];
				$pepNam2 =~ tr/>//d;
				$printedPeps2++;
			}
			else
			{
				$switch = 0;
			}
			$cntr++;
		}
		if($switch == 1)
		{
			print {$allpeps} $line . "\n";
		}
	}
	close PEPS2;
	close $allpeps;

	my $printedPeps = $printedPeps1 + $printedPeps2;

	if($printedPeps1 > 0 && $printedPeps2 > 0)
	{
	
		#system("cat frag*.fa.transdecoder_dir/longest_orfs.pep > allpeps.pep");
		if($printedPeps > 2)
		{
			system("../../../clustalo -i allpeps.pep --distmat-out=peps.dist --percent-id --full --outfmt=vie -o peps.aln --force");
		}
		elsif($printedPeps > 1)
		{
			system("../../../clustalo -i allpeps.pep --full --outfmt=vie -o peps.aln --force");
		}
		else
		{
			next;
		}
	}
	else
	{
		next;
	}
	
	my $f1Count = 0;
	my $f2Count = 0;
	my %maxMatch;
	my %indNams;
	my @lines;
	my $pepCnt = 1;

	if($printedPeps > 2)
	{
		print "Opening peps.dist\n";
		open MATRIX, "peps.dist" or next;
	
		while(<MATRIX>)
		{
			chomp;
			if($. > 1)
			{
				push(@lines, $_);
				my @sps = split(/\s+/);
				$indNams{$pepCnt} = $sps[0];
				my @sps2 = split(/::/, $sps[0]);
				my $ch1 = substr $sps2[1], 4, 1;
			
				if($ch1 == 1)
				{
					$f1Count++;
				}
				if($ch1 == 2)
				{
					$f2Count++;
				}
				$pepCnt++;
			}
		}
	
		close MATRIX;
	}
	else
	{
		$indNams{1} = $pepNam1;
		$indNams{2} = $pepNam2;
		$f1Count = $printedPeps1;
		$f2Count = $printedPeps2;
	}
	#Protlen, matches, mismatches, indels, best-perc
	if($f1Count > 0 && $f2Count > 0)
	{
		$pepCnt = 1;
		my %indSts;
		my %indEds;

		print "opening fr1x.cds\n";
		open FR1, "fr1x.cds";
		while(<FR1>)
		{
			chomp;
			my @sps = split(/\t/);
			$indSts{$pepCnt} = $sps[0];
			$indEds{$pepCnt} = $sps[1];
			$pepCnt++;
		}
		close FR1;
		print "opening fr2x.cds\n";
		open FR2, "fr2x.cds";
		while(<FR2>)
		{
			chomp;
			my @sps = split(/\t/);
			$indSts{$pepCnt} = $sps[0];
			$indEds{$pepCnt} = $sps[1];
			$pepCnt++;
		}
		close FR2;
	
		my %indPairs;
		if($printedPeps > 2)
		{
			my @f1Inds = ();
			my @f2Inds = ();
			my $currInd = 0;

			if($f1Count > 1)
			{
				@f1Inds = (1..$f1Count);
			}
			else
			{
				$f1Inds[0] = 1;
			}
			if($f2Count > 1)
			{
				@f2Inds = (($f1Count + 1)..($f1Count + $f2Count));
			}
			else
			{
				$f1Inds[0] = $f1Count + 1;
			}
			my @f2Inds = (($f1Count + 1)..($f1Count + $f2Count));
			foreach(@lines)
			{
				$currInd++;
				
				my $matchInd = 0;
				my $currMax = 0;

				my @spsL = split(/\s+/, $_);
				my @useInds = @f1Inds;
				if($currInd <= $f1Count)
				{
					@useInds = @f2Inds;
				}
				foreach(@useInds)
				{
					my $tInd = $_;
				
					if($spsL[$tInd] > $currMax)
					{
						$matchInd = $tInd;
						$currMax = $spsL[$tInd];
					}
				}
				$indPairs{$currInd} = $matchInd;
				$currMax = 0;
			}
		}
		else
		{
			$indPairs{1} = 2;
			$indPairs{2} = 1;
		}
		my %peps;
		print "opening allpeps.pep again\n";
		open PEPS, "allpeps.pep" or die $!;
		my $oldID = "";
		my $seq = "";
		
		while(<PEPS>)
		{
			chomp;
			my $line = $_;
			my $subs = substr $_, 0, 1;
			if($subs eq ">")
			{
				if(length($seq) > 2)
				{
					$peps{$oldID} = $seq;
				}
				my @spsA = split(/\s+/);
				$oldID = $spsA[0];
				$oldID =~ tr/>//d;
			}
			else
			{
				$seq = $line;
			}
		}
		$peps{$oldID} = $seq;
		close PEPS;
		
		print "Making directory 'matches'\n";
		system("mkdir matches");
		
		my @keyInds = keys %indPairs;
		my $kcnt = scalar @keyInds;
		my @fInds = (1..$kcnt);
		my @pepFiles = ();
		my @pepAlns = ();

		my @PolyX = ();
		my @PolyY = ();
		
		foreach(@fInds)
		{
			my $ind1 = $_;
			my $ind2 = $indPairs{$ind1};
			
			my $nam1 = $indNams{$ind1};
			my $nam2 = $indNams{$ind2};

			if($ind1 <= $f1Count)
			{
				push(@PolyY, 1);
				push(@PolyY, 1);
				push(@PolyY, 2);
				push(@PolyY, 2);
			}
			else
			{
				push(@PolyY, 2);
				push(@PolyY, 2);
				push(@PolyY, 1);
				push(@PolyY, 1);
			}

			push(@PolyX, $indSts{$ind1});
			push(@PolyX, $indEds{$ind1});
			push(@PolyX, $indSts{$ind2});
			push(@PolyX, $indEds{$ind2});
			
			push(@pepFiles, "matches/" . $ind1 . "_" . $ind2 . ".pep");
			push(@pepAlns, "matches/" . $ind1 . "_" . $ind2 . ".vie");
			
			open my $fout, '>', "matches/" . $ind1 . "_" . $ind2 . ".pep";
			
			print {$fout} ">" . $nam1 . "\n";
			print {$fout} $peps{$nam1} . "\n";
			print {$fout} ">" . $nam2 . "\n";
			print {$fout} $peps{$nam2} . "\n";
			
			close $fout;
		}
		
		my $cntX = 0;

		my @innerCnt = ();
		my @matchCnt = ();
		my @misMatches = ();
		my @totIndel = ();
		my @percs = ();
		my @slens = ();

		my @SubsX0 = ();
		my @SubsY0 = ();
		my @SubsX1 = ();
		my @SubsY1 = ();
	
		foreach(@pepFiles)
		{
			my $sq1Cntr = 0;
			my $sq2Cntr = 0;

			my $file = $_;
			my $outF = $pepAlns[$cntX];
			system("../../../clustalo -i " . $file . " --outfmt=vie -o " . $outF . " --force");
			my $seq1 = "";
			my $seq2 = "";

			my $ind1 = $cntX + 1;
			my $ind2 = $indPairs{$ind1};

			open VIE, $outF or die $!;
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

			$innerCnt[$cntX] = 0;
			$matchCnt[$cntX] = 0;
			$misMatches[$cntX] = 0;
			$totIndel[$cntX] = 0;
			$percs[$cntX] = 0;


			my $len = length($seq1);
			my $gapAccu = 0;

			for(my $i = 0; $i < $len; $i++)
			{
				my $ch1 = substr $seq1, $i, 1;
				my $ch2 = substr $seq2, $i, 1;

				if($ch1 ne "-")
				{
					$sq1Cntr++;
				}
				if($ch2 ne "-")
				{
					$sq2Cntr++;
				}

				if($ch1 eq $ch2)
				{
					if($gapAccu > 0)
					{
						$totIndel[$cntX] += $gapAccu;
						$innerCnt[$cntX] += $gapAccu;
			
						$gapAccu = 0;
					}

					$matchCnt[$cntX]++;
			        $innerCnt[$cntX]++;
				}
				else
				{
					if($ch1 ne "-" && $ch2 ne "-")
		 	        {
						if($gapAccu > 0)
		            	{
	                    	$totIndel[$cntX] += $gapAccu;
	                    	$innerCnt[$cntX] += $gapAccu;

							$gapAccu = 0;
		            	}
						$misMatches[$cntX]++;
	                    $innerCnt[$cntX]++;

						my $s1x = $indSts{$ind1} + ($sq1Cntr * 3);
						my $s2x = $indSts{$ind2} + ($sq2Cntr * 3);
						push(@SubsX0, $s1x);
						push(@SubsX1, $s2x);
						if($ind1 <= $f1Count)
						{
							push(@SubsY0, 1);
							push(@SubsY1, 2);
						}
						else
						{
							push(@SubsY0, 2);
							push(@SubsY1, 1);
						}
					}
					else
					{
						$gapAccu++;
					}
				}
			}		
			$percs[$cntX] = ($matchCnt[$cntX] / ($misMatches[$cntX] + $matchCnt[$cntX])) * 100;
			$slens[$cntX] = $sq1Cntr;
			$cntX++;
		}
		open my $statsOut, '>', "prot.stats" or die $!;
	
		#print {$statsOut} "Protlen\tmatches\tmismatches\tindels\tbest-perc\n";

		for(my $f = 0; $f < $cntX; $f++)
		{
			print {$statsOut} $slens[$f] . "\t" . $matchCnt[$f] . "\t" . $misMatches[$f] . "\t" . $totIndel[$f] . "\t" . $percs[$f] . "\n";
		}

		close $statsOut;

		open my $polyOut, '>', "prot.polys" or die $!;

		my $polyLen = scalar @PolyY;

		for(my $f = 0; $f < $polyLen; $f++)
		{
			print {$polyOut} $PolyX[$f] . "\t" . $PolyY[$f] . "\n";
		}

		close $polyOut;

		open my $linkOut, '>', "prot.subs" or die $!;

		my $linkLen = scalar @SubsX0;

		print {$linkOut} "x0\ty0\tx1\ty1\n";

		for(my $f = 0; $f < $linkLen; $f++)
		{
			print {$linkOut} $SubsX0[$f] . "\t" . $SubsY0[$f] . "\t" . $SubsX1[$f] . "\t" . $SubsY1[$f] . "\n";
		}

		close $linkOut;	
	}

	chdir("../../../");
	$success = 1;
}
