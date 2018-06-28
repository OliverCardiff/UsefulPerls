use strict;
use warnings;

my $READ_DEPTH = 10;

my $skips = 0;

my $seq = "";

my $scaff = 0;

while(<>)
{
	chomp;
	my @sps = split(/\t/);
	if(scalar(@sps) > 4)
	{
	
		if($skips == 0)
		{

			my $mess = $sps[4];

			my $chosen = "";

			if($sps[3] < $READ_DEPTH)
			{
				$chosen = "N";
			}
			else
			{
				$sps[4] =~ tr/atgc/TACG/d;

				my @letters = qw(A T G C);
				my @freqs = qw(0 0 0 0);
				my $mxFrq = 0;
				$chosen = $sps[2];
		
				for(my $it = 0; $it < 4; $it++)
				{
					my $l = $letters[$it];
					$freqs[$it] = () = $sps[4] =~ /$l/g;
					if($freqs[$it] > $mxFrq)
					{
						$mxFrq = $freqs[$it];
						$chosen = $l;
					}
				}
			
				@letters = qw(- \+);
				@freqs = qw(0 0);
				$mxFrq = 0;
				my $sign = "";
		
				for(my $it = 0; $it < 2; $it++)
				{
					my $l = $letters[$it];
					$freqs[$it] = () = $sps[4] =~ /$l/g;
					if($freqs[$it] > $mxFrq)
					{
						$mxFrq = $freqs[$it];
						$sign = $l;
					}
				}
			
				if($sign ne "" && $mxFrq > 2)
				{
					@letters = qw(1 2 3 4 5 6 7 8 9);
					@freqs = qw(0 0 0 0 0 0 0 0 0);
					$mxFrq = 0;
					my $size = 1;
		
					for(my $it = 0; $it < 9; $it++)
					{
						my $l = $letters[$it];
						$freqs[$it] = () = $sps[4] =~ /$l/g;
						if($freqs[$it] > $mxFrq)
						{
							$mxFrq = $freqs[$it];
							$size = $l;
						}
					}
				
					if($sign eq "\\+")
					{
						my $find = $sign . $size;

						my ($str) = ($sps[4] =~ /$find[0-9]*(.*)\./);

						$chosen = substr $str, 0, $size;
						#print "In the plus sign : " . $chosen . "\n";
					
					}
					else
					{
						$skips += $size - 1;
						$chosen = "";
					}
				}
			}

			$seq .= $chosen;

			#Sanity checker
			#if($chosen ne $sps[2])
			#{
			#	print $sps[2] . "\t" . $chosen . "\n";
			#}
		}
		else
		{
			$skips--;
		}
	}
	else
	{
		if($seq ne "")
		{
			printseqs();
		}
	}
}

if($seq ne "")
{
	printseqs();
}

sub printseqs {

	print ">My_consensus_plz_rename_me_" . $scaff . "\n";

	for(my $i = 0; $i < length($seq) - 100; $i += 100)
	{
		print substr($seq, $i, 100) . "\n";
	}

	my $ln = length($seq) % 100;

	print substr($seq, length($seq) - $ln, $ln) . "\n";

	$seq = "";
	$scaff++;
}
