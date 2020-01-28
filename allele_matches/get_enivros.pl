#usage perl get_enviros.pl stub_prostats.txt stub2_prostats.txt  stub.tblout stub2.tblout select_pfams.txt
use strict;
use warnings;

my %pepLines1;
my %pepLines2;

my %pepPfam;
my %pepPfam2;

my @pfams;
my %pfamCnt;
my %pfamCnt2;

my %pfamPeps1;
my %pfamPeps2;

my $headerLine = "";

open PROSTAT, $ARGV[0] or die $!;

while(<PROSTAT>)
{
	chomp;
	if($. == 1)
	{
		$headerLine = $_;
	}
	my @sps = split(/\t/);
	$pepLines1{$sps[0]} = $_;
}

close PROSTAT;

open PROSTAT2, $ARGV[1] or die $!;

while(<PROSTAT2>)
{
        chomp;
        my @sps = split(/\t/);
        $pepLines2{$sps[0]} = $_;
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

			push(@{$pfamPeps1{$pfid[0]}}, $sps[2]);

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

						push(@{$pfamPeps2{$pfid[0]}}, $sps[2]);
                }
                $oldID = $sps[2];
        }
}

close TBLOUT2;

open TARGET, $ARGV[4] or die $!;

my @chosenOnes;

while(<TARGET>)
{
	chomp;
	my @sps = split(/\t/);

	push(@chosenOnes, $sps[0]);	
}

close TARGET;

foreach(@chosenOnes)
{
	my $pf = $_;
	my @pfpep1 = @{$pfamPeps1{$pf}};
	my @pfpep2 = @{$pfamPeps2{$pf}};
	print "Writing enviro/ling_" . $pf . ".txt\n";

	open my $fout1, '>', "enviro/ling_" . $pf . ".txt" or die $!;
	print {$fout1} $headerLine . "\n";	
	foreach(@pfpep1)
	{
		my $pp = $_;
		print {$fout1} $pepLines1{$pp} . "\n";
	}
	close $fout1;

	print "Writing enviro/rube_" . $pf . ".txt\n";

	open my $fout2, '>', "enviro/rube_" . $pf . ".txt" or die $!;
	print {$fout2} $headerLine . "\n";	
	foreach(@pfpep2)
	{
		my $pp = $_;
		print {$fout2} $pepLines2{$pp} . "\n";
	}
	close $fout2;
}

