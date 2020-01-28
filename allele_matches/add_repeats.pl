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
	system("mkdir repeats");
	system("blastn -query frag1.fa -db ../../../" . $stub . "_repeats/$stub.fa.classified -outfmt 6 -evalue 5e-05 > repeats/fr1.blast");
	system("blastn -query frag2.fa -db ../../../" . $stub . "_repeats/$stub.fa.classified -outfmt 6 -evalue 5e-05 > repeats/fr2.blast");

	open my $rep1, '>', "fr1.rep" or next;
	open REP1, "repeats/fr1.blast" or next;

	while(<REP1>)
	{
		my @sps = split(/\t/);
		
		if($sps[3] > 80)
		{
			if($sps[8] > $sps[9])
			{
				my $temp = $sps[9];
				$sps[9] = $sps[8];
				$sps[8] = $temp;
			}
			print {$rep1} $sps[8] . "\t". $sps[9] . "\n"
		}
	}

	close REP1;
	close $rep1;

	open my $rep2,  '>', "fr2.rep" or next;
	open REP2, "repeats/fr2.blast" or next;

	while(<REP2>)
	{
		my @sps = split(/\t/);
		
		if($sps[3] > 80)
		{
			if($sps[8] > $sps[9])
			{
				my $temp = $sps[9];
				$sps[9] = $sps[8];
				$sps[8] = $temp;
			}
			print {$rep2} $sps[8] . "\t". $sps[9] . "\n";
		}
	}

	close REP2;
	close $rep2;

	$success = 1;
	chdir("../../../");
}
