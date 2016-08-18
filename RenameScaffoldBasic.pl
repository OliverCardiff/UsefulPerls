use strict;
use warnings;

my %AG_GHash;

open AG_G, $ARGV[0] or die $!;
open my $fout, '>', $ARGV[1] or die $!;
open my $lift, '>', $ARGV[2] or die $!;

my $cnt = 1;
while(<AG_G>)
{
	chomp;
	my $line = $_;
	
	my $subs = substr $_, 0, 1;
	
	if($subs eq ">")
	{
		my $line2 = $line;
		$line = ">scaffold_" . $cnt;
		print {$lift} $line . "\t". $line2 . "\n";
		$cnt++;
	}
	
	print {$fout} $line . "\n";
}

close $fout;
close AG_G;
close $lift;
