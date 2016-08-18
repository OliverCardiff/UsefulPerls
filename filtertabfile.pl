use strict;
use warnings;

my $indexPos = $ARGV[1];

my %scaffHash;

open LIST, $ARGV[0] or die $!;

while(<LIST>)
{
	chomp;
	my $scaff = $_;
	$scaffHash{$scaff} = 1;
}

close LIST;

open TABFILE, $ARGV[2] or die $!;
open my $filt, '>', $ARGV[3] or die $!;

my $switch = 0;

while(<TABFILE>)
{
	chomp;
	my $line = $_;
	my @sps = split(/\s+/);
	
	if(exists $scaffHash{$sps[$indexPos]} || exists $scaffHash{$sps[$indexPos + 1]})
	{
		$switch = 0;
	}
	else
	{
		$switch = 1;
	}
	
	if($switch == 1)
	{
		print {$filt} $line  . "\n";
	}
}

close $filt;
close TABFILE;
