#usage perl make_go_mats.pl stub_prostats.txt stub2_prostats.txt  stub.tblout stub2.tblout [out]pfam_counts
use strict;
use warnings;

my %pepPerc;
my %pepPerc2;

my %pepPfam;
my %pepPfam2;

my @pfams;
my %pfamCnt;
my %pfamCnt2;

my %pfamAccu1;
my %pfamAccu2;

open PROSTAT, $ARGV[0] or die $!;

while(<PROSTAT>)
{
	chomp;
	my @sps = split(/\t/);
	$pepPerc{$sps[0]} = $sps[7];
}

close PROSTAT;

open PROSTAT2, $ARGV[1] or die $!;

while(<PROSTAT2>)
{
        chomp;
        my @sps = split(/\t/);
        $pepPerc2{$sps[0]} = $sps[7];
}

close PROSTAT2;

open TBLOUT, $ARGV[2] or die $!;
my $oldID = "";
while(<TBLOUT>)
{
	chomp;
	my $subs = substr $_,0,1;
	if($. > 3 && $subs ne "#")
	{
		my @sps = split(/\s+/);
		if($oldID ne $sps[2])
		{
			my @pfid = split(/\./, $sps[1]);
			$pepPfam{$sps[2]} = $pfid[0];
			if(!exists $pfamCnt{$pfid[0]})
			{
				$pfamCnt{$pfid[0]} = 1;
				push(@pfams, $pfid[0]);
			}
			else
			{
				$pfamCnt{$pfid[0]}++;
			}
			if(!exists $pfamAccu1{$pfid[0]})
			{
				$pfamAccu1{$pfid[0]} = $pepPerc{$sps[2]};
			}
			else
			{
				$pfamAccu1{$pfid[0]} += $pepPerc{$sps[2]};
			}
		}
		$oldID = $sps[2];
	}
}

close TBLOUT;

open TBLOUT2, $ARGV[3] or die $!;
$oldID = "";
while(<TBLOUT2>)
{
        chomp;
        my $subs = substr $_,0,1;
        if($. > 3 && $subs ne "#")
        {
                my @sps = split(/\s+/);
                if($oldID ne $sps[2])
                {
                        my @pfid = split(/\./, $sps[1]);
                        $pepPfam2{$sps[2]} = $pfid[0];
                        if(!exists $pfamCnt2{$pfid[0]})
                        {
                                $pfamCnt2{$pfid[0]} = 1;
                        }
                        else
                        {
                                $pfamCnt2{$pfid[0]}++;
                        }
			if(!exists $pfamAccu2{$pfid[0]})
                        {
                                $pfamAccu2{$pfid[0]} = $pepPerc2{$sps[2]};
                        }
                        else
                        {
                                $pfamAccu2{$pfid[0]} += $pepPerc2{$sps[2]};
                        }
                }
                $oldID = $sps[2];
        }
}

close TBLOUT2;

open my $fout_ref, '>', $ARGV[4] or die $!;

print {$fout_ref} "pfamID\tCount1\tCount2\tCountSum\tAvgDiverge1\tAvgDiverge2\tAvgDivBoth\n";

foreach(@pfams)
{
	my $pf = $_;
	
	my $cnt1 = 0;
	my $cnt2 = 0;
	my $p1 = 100;
	my $p2 = 100;
	if(exists $pfamCnt{$pf})
	{
		$cnt1 = $pfamCnt{$pf};
		$p1 = $pfamAccu1{$pf}/$cnt1;
	}
	if(exists $pfamCnt2{$pf})
        {
                $cnt2 = $pfamCnt2{$pf};
		$p2 = $pfamAccu2{$pf}/$cnt2;
        }
	if($cnt2 > 1 && $cnt1 > 1)
	{
		my $tot = $cnt1 + $cnt2;
		my $avg = ($p2 +$p1)/2;
		print {$fout_ref} "$pf\t$cnt1\t$cnt2\t$tot\t$p1\t$p2\t$avg\n";
	}
}

close $fout_ref;
