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

#print join("ppp", @) . "\n";

my $spanID = 1;

foreach(@spansL)
{
	my $span = $_;
	my $line = $binHash{$span};
	my @sps = split(/\t/, $line);
	my $len = scalar @sps;
	my @binaries = @sps[3..($len-1)];

	my $basicStr = $spanID . "\t" . $span . "\t" . join("\t", @binaries);
	$spanID++;

	print $basicStr . "\n";
}


