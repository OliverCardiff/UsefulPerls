use strict;
use warnings;

my %scaffs;

open LIST, $ARGV[0] or die $!;

while(<LIST>)
{
	chomp;

	if(exists $scaffs{$_})
	{
		$scaffs{$_}++;
	}
	else
	{
		$scaffs{$_} = 1;
	}
}

close LIST;

open IDX, $ARGV[1] or die $!;

while(<IDX>)
{
	chomp;
	my $line = $_;
	my @sps = split(/\t/); 

	if(exists $scaffs{$sps[0]})
	{
		if($scaffs{$sps[0]} == 3)
		{
			print $line . "\n";
		}
	}
}

close IDX;

