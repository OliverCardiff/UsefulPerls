#usage perl get_nets.pl namvec_full.txt namvec_single.txt namvec_sparse.txt peps.blast output.txt

use strict;
use warnings;

my %scaffLens;
my %scaffSpans;
my %scaffSelfs;
my %scaffPeps;
my %spans;

my %spanInternal;
my %spanExternal;
my %spanSpans;

my $pepCount = 0;

for(my $i = 0; $i < 3; $i++)
{
	my $file = $ARGV[$i];
	
	open SPAN, $file or die $!;
	
	while(<SPAN>)
	{
		if($. > 1)
		{
			my @sps = split(/\t/);
			$scaffSpans{$sps[0]} = $sps[2];
			$spans{$sps[2]} = $sps[1];
		}
	}
	
	close SPAN;
}

open PEPS, $ARGV[3] or die $!;

while(<PEPS>)
{
	chomp;
	my @sps = split(/\t/);
	my @scsps = split(/::/, $sps[0]);
	my @scspsB = split(/::/, $sps[1]);
	my @sc2 = split(/\|/, $scsps[1]);
	#print $scspsB[1] . "\t" . $scsps[1] . "\n";
	my $len = substr $sc2[1], 4;
	$pepCount++;
	
	if(exists $scaffSpans{$scsps[1]} && exists $scaffSpans{$scspsB[1]})
	{
		my $spanA = $scaffSpans{$scsps[1]};
		my $spanB = $scaffSpans{$scspsB[1]};
		
		if(exists $spanSpans{$spanA}{$spanB})
		{
			$spanSpans{$spanA}{$spanB}++;
		}
		else
		{
			$spanSpans{$spanA}{$spanB} = 1;
		}
	}
	$scaffLens{$scsps[1]} = $len;
	if(exists $scaffPeps{$scsps[1]})
	{
		$scaffPeps{$scsps[1]}++;
	}
	else 
	{
		$scaffPeps{$scsps[1]} = 1;
	}
	
	if($scspsB[1] eq $scsps[1])
	{
		if(exists($scaffSelfs{$scsps[1]}))
		{
			$scaffSelfs{$scsps[1]}++;
		}
		else
		{
			$scaffSelfs{$scsps[1]} = 1;
		} 
	}
	else
	{
		if(defined $scaffSpans{$scsps[1]})
		{
			my $spanA = $scaffSpans{$scsps[1]};
			if(exists $spanExternal{$spanA})
			{
				$spanExternal{$spanA}++;
			}
			else
			{
				$spanExternal{$spanA} = 1;
			}
		}
	}
}

close PEPS;

my $accuIn = 0;
my $accuOut = 0;
my @allScaffs = keys %scaffLens;

foreach(@allScaffs)
{
	my $scaff = $_;
	my $peps = $scaffPeps{$scaff};
	my $ins = $scaffSelfs{$scaff};
	if(!defined $ins)
	{
		$ins = 0;
	}
	my $outs = $peps - $ins;
	$accuIn += $ins;
	$accuOut += $outs;
}

my $avgRatio = $accuIn / ($accuIn + $accuOut);

my @spanAll = keys %spans;

foreach(@spanAll)
{
	my $spanA = $_;
	print $spanA . "\t" . $spanA . "\n";
}

=retired
my $spanAccu = 0;
my $rateAccu = 0;

foreach(@spanAll)
{
	$spanAccu++;

	my $span = $_;
	my $in = $spanSpans{$span}{$span};
	my $out = $spanExternal{$span};
	my $rate = $in / ($in + $out);
	my $factor = $rate / $avgRatio;

	$rateAccu += $factor;
	
	print $span . "\t" . $in . "\t" . $out . "\t" . $rate . "\t" . $avgRatio . "\t" . $factor . "\n";
}

my $finRate = $rateAccu / $spanAccu;
print "\nFinal Rate: " . $finRate . "\n";
=cut