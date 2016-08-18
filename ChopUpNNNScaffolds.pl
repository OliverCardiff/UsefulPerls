use strict;
use warnings;

my %ch12Hash;

open CH12, $ARGV[0] or die $!;
open my $fout, '>', $ARGV[1] or die $!;

my $nextID = "";
my $fragID = 1;
my $Nflag = 0;

while(<CH12>)
{
	chomp;
	my $line = $_;
	my $subs = substr $_, 0, 1;
	
	if($subs eq ">")
	{
		$nextID = substr $line, 0, 26;
		$fragID = 1;
		print {$fout} $nextID . "-" . $fragID . "\n";
		$fragID++;
	}
	else
	{
		if(index($_, "NNNNNNNNNNNNNNNNNNNNNNNNNN") != -1)
		{
			if($Nflag == 0)
			{
				print {$fout} $nextID . "-" . $fragID . "\n";
				$fragID++;
				$Nflag = 1;
			}
		}
		else
		{
			$Nflag = 0;
			
			print {$fout} $line . "\n";
		}
	}	
}

close $fout;
close CH12;
