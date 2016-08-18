use strict;
use warnings;

my %KO_1;
my %KO_2;
my %KO_3;
my %KO_4;

open FOAM, $ARGV[0] or die $!;

while(<FOAM>)
{
	if($. != 1)
	{
		chomp;
		my @sps = split(/\t/);
		$KO_1{$sps[4]} = $sps[0];
		$KO_2{$sps[4]} = $sps[1];
		$KO_3{$sps[4]} = $sps[2];
		$KO_4{$sps[4]} = $sps[3];
	}
}

close FOAM;

print "Read Foam table\n";
my %scaffKO;
my $oldId = "";


open HMM, $ARGV[1] or die $!;

while(<HMM>)
{
	chomp;
	my $subs = substr $_, 0, 1;
	if($subs ne "#")
	{
		my @sps = split(/\s+/);
		my @sps2 = split(/::/, $sps[2]);
		my $ko = substr $sps[0], 3, 6;

		my $scaff = $sps2[1];
		
		if(!($oldId eq $scaff))
		{
			$scaffKO{$scaff} = $ko;
		}
		
		$oldId = $scaff;
	}
	
}

close HMM;
print "Read HMM table\n";
my %scaffCnt;

open IDX, $ARGV[2] or die $!;

while(<IDX>)
{
	chomp;
	my @sps = split(/\t/);
	
	if(exists $scaffKO{$sps[0]})
	{
		$scaffCnt{$sps[0]} = $sps[2]/$sps[1];
	}
}

close IDX;
print "Read IDX stats\n";

open my $lvl1, '>', $ARGV[3] . "_lvl1.txt" or die $!;
open my $lvl2, '>', $ARGV[3] . "_lvl2.txt" or die $!;
open my $lvl3, '>', $ARGV[3] . "_lvl3.txt" or die $!;
open my $lvl4, '>', $ARGV[3] . "_lvl4.txt" or die $!;

my %L_Hash;

my @keys = keys %scaffKO;

foreach(@keys)
{
	my $sc = $_;
	my $ko = $scaffKO{$sc};
	my $l1 = $KO_1{$ko};
	my $l2 = $KO_2{$ko};
	my $l3 = $KO_3{$ko};
	my $l4 = $KO_4{$ko};
	
	if(defined $l1 && !($l1 eq ""))
	{
		if(exists $L_Hash{"1"}{$l1}{"reads"})
		{
			$L_Hash{"1"}{$l1}{"reads"} += $scaffCnt{$sc};
		}
		else
		{
			$L_Hash{"1"}{$l1}{"reads"} = $scaffCnt{$sc};
		}
		
		if(exists $L_Hash{"1"}{$l1}{"genes"})
		{
			$L_Hash{"1"}{$l1}{"genes"} += 1;
		}
		else
		{
			$L_Hash{"1"}{$l1}{"genes"} = 1;
		}
	}
	if(defined $l2 && !($l2 eq ""))
	{
		if(exists $L_Hash{"2"}{$l2}{"reads"})
		{
			$L_Hash{"2"}{$l2}{"reads"} += $scaffCnt{$sc};
		}
		else
		{
			$L_Hash{"2"}{$l2}{"reads"} = $scaffCnt{$sc};
		}
		
		if(exists $L_Hash{"2"}{$l2}{"genes"})
		{
			$L_Hash{"2"}{$l2}{"genes"} += 1;
		}
		else
		{
			$L_Hash{"2"}{$l2}{"genes"} = 1;
		}
	}
	if(defined $l3 && !($l3 eq ""))
	{
		if(exists $L_Hash{"3"}{$l3}{"reads"})
		{
			$L_Hash{"3"}{$l3}{"reads"} += $scaffCnt{$sc};
		}
		else
		{
			$L_Hash{"3"}{$l3}{"reads"} = $scaffCnt{$sc};
		}
		
		if(exists $L_Hash{"3"}{$l3}{"genes"})
		{
			$L_Hash{"3"}{$l3}{"genes"} += 1;
		}
		else
		{
			$L_Hash{"3"}{$l3}{"genes"} = 1;
		}
	}
	if(defined $l4 && !($l4 eq ""))
	{
		if(exists $L_Hash{"4"}{$l4}{"reads"})
		{
			$L_Hash{"4"}{$l4}{"reads"} += $scaffCnt{$sc};
		}
		else
		{
			$L_Hash{"4"}{$l4}{"reads"} = $scaffCnt{$sc};
		}
		
		if(exists $L_Hash{"4"}{$l4}{"genes"})
		{
			$L_Hash{"4"}{$l4}{"genes"} += 1;
		}
		else
		{
			$L_Hash{"4"}{$l4}{"genes"} = 1;
		}
	}
}
print "Built Hash of all data\n";

my @l1_keys = keys %{$L_Hash{"1"}};
my @l2_keys = keys %{$L_Hash{"2"}};
my @l3_keys = keys %{$L_Hash{"3"}};
my @l4_keys = keys %{$L_Hash{"4"}};

foreach(@l1_keys)
{
	my $nam = $_;
	print {$lvl1} $nam . "\t" . $L_Hash{"1"}{$nam}{"genes"} . "\t" . $L_Hash{"1"}{$nam}{"reads"} . "\n";
}
print "Wrote Tier 1 funcs to: " . $ARGV[3] . "_lvl1.txt\n";
foreach(@l2_keys)
{
	my $nam = $_;
	print {$lvl2} $nam . "\t" . $L_Hash{"2"}{$nam}{"genes"} . "\t" . $L_Hash{"2"}{$nam}{"reads"} . "\n";
}
print "Wrote Tier 2 funcs to: " . $ARGV[3] . "_lvl2.txt\n";
foreach(@l3_keys)
{
	my $nam = $_;
	print {$lvl3} $nam . "\t" . $L_Hash{"3"}{$nam}{"genes"} . "\t" . $L_Hash{"3"}{$nam}{"reads"} . "\n";
}
print "Wrote Tier 3 funcs to: " . $ARGV[3] . "_lvl3.txt\n";
foreach(@l4_keys)
{
	my $nam = $_;
	print {$lvl4} $nam . "\t" . $L_Hash{"4"}{$nam}{"genes"} . "\t" . $L_Hash{"4"}{$nam}{"reads"} . "\n";
}
print "Wrote Tier 4 funcs to: " . $ARGV[3] . "_lvl4.txt\n";
close $lvl1; close $lvl2; close $lvl3; close $lvl4;
