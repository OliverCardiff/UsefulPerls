#usage perl pairs_filt_and_blast.pl <final.mat> <genome.fa> <out.prefix> <integer - smallest match>
use strict;
use warnings;

my %pairs;

my $smallest = $ARGV[3];
my $stub = $ARGV[2];

open FINAL, $ARGV[0] or die $!;

while(<FINAL>)
{
	chomp;
	if($. > 1)
	{

		my @sps = split(/\t/);

		my $smol  = $sps[2];

		if($sps[3] < $smol)
		{
			$smol = $sps[3];
		}

		if($smol > $smallest && !(exists $pairs{$sps[1]}))
		{
			$pairs{$sps[0]}{$sps[1]} = 1;
		}
	}
}


close FINAL;

my @primaries = keys %pairs;
my $scCnt = scalar @primaries;

print "Found $scCnt pair matches above $smallest length!\n";

my %genome;

open GENOME, $ARGV[1] or die $!;

my $oldID = "";
my $seq = "";

while(<GENOME>)
{
	chomp;
	my $line = $_;
	my $subs = substr $_, 0, 1;

	if($subs eq ">")
	{
		if(length($seq) > 2)
		{
			$genome{$oldID} = $seq;
		}
		$oldID = $line;
		$oldID =~ tr/>//d;
		$seq = "";
	}
	else
	{
		$seq .= $line;
	}
}

close GENOME;

my $cnt = 1;

foreach(@primaries)
{
	my $sc1 = $_;
	my @sc_tr = split(/\|/, $sc1);
	my $sc_file = $sc_tr[0];

	system("mkdir $stub". "_" . $cnt);

	open my $fout1, '>', $stub . "_$cnt/$sc_file.fa" or die $!;

	print {$fout1} ">$sc1\n";
	print {$fout1} $genome{$sc1} . "\n";
	
	close $fout1;

	open my $all, '>', $stub . "_$cnt/all.fa" or die $!;
	print {$all} ">$sc1\n";
	print {$all} $genome{$sc1} . "\n";

	my @matches = keys %{$pairs{$sc1}};

	foreach(@matches)
	{
		my $sc2 = $_;
		my @sc_tr2 = split(/\|/, $sc2);
		my $sc2_file = $sc_tr2[0];

		open my $fout2, '>', $stub . "_$cnt/$sc2_file.fa" or die $!;

		print {$fout2} ">$sc2\n";
		print {$fout2} $genome{$sc2} . "\n";

		close $fout2;

		print {$all} ">$sc2\n";
		print {$all} $genome{$sc2} . "\n"
	}

	close $all;

	system("makeblastdb -in $stub" . "_$cnt/all.fa -dbtype nucl -parse_seqids");
	system("blastn -query $stub" . "_$cnt/all.fa -db $stub" . "_$cnt/all.fa -outfmt 5 -out $stub" . "_$cnt/$sc_file.blast.xml");
	$cnt++;	
}

