#usage perl make_over_reps.pl namvec.txt pfamgo.txt grac.tblout pfam_sheet.out go_sheet.out

use strict;
use warnings;

my %scafSpans;
my %binHash;
my %pfamAccu1;
my %goAccu1;
my %ABids;
my $oldID = "";
my @spansL;

my @inds = qw(CWMC_1-1	CWMC_1-2	CWMC_1-3	CWMC_1-4	CWMC_1-5	CWMC_1-6	CWMC_1-7	CWMC_1-8	CWMC_2-10	CWMC_2-1	CWMC_2-2	CWMC_2-3	CWMC_2-4	CWMC_2-5	CWMC_2-6	CWMC_2-7	CWMC_2-8	CWMC_2-9	CWMC_3-10	CWMC_3-2	CWMC_3-3	CWMC_3-4	CWMC_3-5	CWMC_3-6	CWMC_3-7	CWMC_3-8	CWMC_3-9	CWMM_1-10	CWMM_1-1	CWMM_1-2	CWMM_1-3	CWMM_1-4	CWMM_1-5	CWMM_1-6	CWMM_1-7	CWMM_1-8	CWMM_1-9	CWMM_2-10	CWMM_2-1	CWMM_2-2	CWMM_2-3	CWMM_2-4	CWMM_2-5	CWMM_2-6	CWMM_2-7	CWMM_2-8	CWMM_2-9	CWMM_3-10	CWMM_3-1	CWMM_3-2	CWMM_3-3	CWMM_3-4	CWMM_3-5	CWMM_3-6	CWMM_3-7	CWMM_3-8	CWMM_3-9);

open AB, "AB_ids.txt" or die $!;

while(<AB>)
{
	chomp;
	my @sps = split(/\s/);
	$ABids{$sps[0]} = $sps[1];
}

close AB;

foreach(@ARGV)
{
	my $file = $_;
	open NAMVEC, $file or die $!;

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
				push(@spansL, $sps[2]);
			}

			$scafSpans{$sps[0]} = $sps[2];
			#print $sps[0] . "\n";
			$oldID = $sps[2];
		}
	}

	close NAMVEC;
}

my @ablist;
my @namlist;

foreach(@inds)
{
	my $ind = $_;
	if(exists $ABids{$ind})
	{
		push(@namlist, $ind);
		push(@ablist, $ABids{$ind});
	}
}

#print join("ppp", @) . "\n";

my $spanID = 1;
print "SpanID	SpanLength	WormCount	" . join("\t", @namlist) ."\n";
print "x\tx\tx\t" . join("\t", @ablist) . "\n";

foreach(@spansL)
{
	my $span = $_;
	my $line = $binHash{$span};
	my @sps = split(/\t/, $line);
	my @binaries = @sps[3..43];

	my $basicStr = $spanID . "\t" . $span . "\t" . join("\t", @binaries);
	$spanID++;

	print $basicStr . "\n";
}


