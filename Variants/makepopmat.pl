use strict;
use warnings;

my %scaffRoll;

open ROLL, $ARGV[0] or die $!;

while(<ROLL>)
{
	chomp;
	
	my @sps = split(/\t/);
	
	push(@{$scaffRoll{$sps[0]}}, $sps[2]);
}

close ROLL;

my %indHash;
my %fullHash;

open MAT, $ARGV[1] or die $!;

while(<MAT>)
{
	chomp;
	my @sps = split(/\t/);
	
	$indHash{$sps[1]} = 1;
	
	$fullHash{$sps[0]}{$sps[2]}{$sps[1]} = 1;
}

close MAT;

my @inds;

my $cnt = 1;
while($cnt < 46)
{
	push (@inds, $cnt);
	$cnt++;
}

open my $fout, '>', $ARGV[2] or die $!;

print {$fout} "Scaffold\tlocus\t" . join("\t", @inds) . "\tRMean\n";

my @scfs = keys %fullHash;

foreach(@scfs)
{
	my $scaff = $_;
	
	my @locs = keys %{$fullHash{$scaff}};
	
	foreach(@locs)
	{
		my $loc = $_;
		my @inds2 = @inds;
		print {$fout} $scaff . "\t" . $loc;
		foreach(@inds2)
		{
			my $ind = $_;
			
			if(exists $fullHash{$scaff}{$loc}{$ind})
			{
				print {$fout} "\t1";
			} else {
				print {$fout} "\t0";
			} 
		}
		
		my $fif = int($loc/50);
		
		print {$fout} "\t" . $scaffRoll{$scaff}[$fif] . "\n";
	}
	print "Did $scaff \n";
}

close $fout;
