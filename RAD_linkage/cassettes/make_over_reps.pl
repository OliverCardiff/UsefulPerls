#usage perl make_over_reps.pl namvec.txt pfamgo.txt grac.tblout pfam_sheet.out go_sheet.out

use strict;
use warnings;

my %pf_cnts;
my %go_cnts;
my %spn_pf_cnts;
my %spn_go_cnts;
my $total_pfs = 0;
my $total_gos = 0;

my %spanPfs;
my %scafSpans;
my %pfamGo;
my %binHash;
my %pfamAccu1;
my %goAccu1;

my $oldID = "";

open NAMVEC, $ARGV[0] or die $!;

while(<NAMVEC>)
{
	chomp;
	my $line = $_;
	if($. > 1)
	{
		my @sps = split(/\t/);
		
		if($oldID ne $sps[2])
		{
			$binHash{$sps[2]} = $line;
		}

		$scafSpans{$sps[0]} = $sps[2];
		#print $sps[0] . "\n";
		$oldID = $sps[2];
	}
}

close NAMVEC;

open PFGO, $ARGV[1] or die $!;

while(<PFGO>)
{
        chomp;
        my @sps = split(/\s+/);
        my $len = scalar @sps;
        
        my $goterm = $sps[$len-1];
        my @pfam = split(/:/, $sps[0]);
        push(@{$pfamGo{$pfam[1]}}, $goterm);
}

close PFGO;

open TBLOUT, $ARGV[2] or die $!;
$oldID = "";
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

			my @scsps = split(/\|/, $sps[2]);
			my $scaff = $scsps[0] . "\|" . $scsps[1];
			#print "$scaff \n";
			if(exists $pf_cnts{$pfid[0]})
			{
				$pf_cnts{$pfid[0]}++;
			}
			else
			{
				$pf_cnts{$pfid[0]} = 1;
			}
			my $span = $scafSpans{$scaff};
			$total_pfs++;

			if(defined $span)
			{

				if(exists $spn_pf_cnts{$span})
				{
					$spn_pf_cnts{$span}++;
				}				
				else {
					$spn_pf_cnts{$span} = 1;
				}

				#print "did it \n";
				#print $sps[0] . "\t" . $sps[1] . "\n";
				if(!exists $pfamAccu1{$span}{$pfid[0]})
				{
					$pfamAccu1{$span}{$pfid[0]} = 1;
				}
				else
				{
					$pfamAccu1{$span}{$pfid[0]} += 1;
				}
				if(exists $pfamGo{$pfid[0]})
				{
					my @gos = @{$pfamGo{$pfid[0]}};

					foreach(@gos)
					{
						my $go = $_;
						$total_gos++;

						if(exists $go_cnts{$go})
						{
							$go_cnts{$go}++;
						}
						else
						{
							$go_cnts{$go} = 1;
						}
						
						if(exists $spn_go_cnts{$span})
						{
							$spn_go_cnts{$span}++;
						}				
						else {
							$spn_go_cnts{$span} = 1;
						}

						if(!exists $goAccu1{$span}{$go})
						{
							$goAccu1{$span}{$go} = 1;
						}
						else
						{
							$goAccu1{$span}{$go} += 1;
						}
					}
				}
			}
		}
		$oldID = $sps[2];
	}
}

close TBLOUT;

open my $foutPF, '>', $ARGV[3] or die $!;
open my $foutGO, '>', $ARGV[4] or die $!;

my @spansL = keys %binHash;
my $spanID = 1;
foreach(@spansL)
{
	my $span = $_;
	my $line = $binHash{$span};
	my @sps = split(/\t/, $line);
	my @binaries = @sps[3..43];

	my @pfams = keys %{$pfamAccu1{$span}};
	my @gos = keys %{$goAccu1{$span}};

	my $basicStr = $spanID . "\t" . $span . "\t" . join("\t", @binaries);

	foreach(@pfams)
	{
		my $pfam = $_;
		my $pfCnt = $pfamAccu1{$span}{$pfam};

		my $R1_spnPF = $pfCnt/$spn_pf_cnts{$span};
		
		my $R2_pfSet = $pf_cnts{$pfam} / $total_pfs;

		my $B_Fact = $R1_spnPF / $R2_pfSet;
		

		if(($B_Fact > 15 && $pfCnt > 1) || ($pfCnt > 4))
		{
			print {$foutPF} $basicStr . "\t" . $pfam . "\t" . $pfCnt . "\t" . $B_Fact . "\n";
		}
	}
	foreach(@gos)
	{
		my $go = $_;
		my $goCnt = $goAccu1{$span}{$go};

		my $R1_spnGO = $goCnt/$spn_go_cnts{$span};
		
		my $R2_goSet = $go_cnts{$go} / $total_gos;

		my $B_Fact = $R1_spnGO / $R2_goSet;

		#if($goCnt > 1)
		#{
			print {$foutGO} $basicStr . "\t" . $go . "\t" . $goCnt . "\t" . $B_Fact ."\n";
		#}
	}
}

close $foutPF;
close $foutGO;
