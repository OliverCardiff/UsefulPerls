# 0 Fasta.fai 1 Genome 2 VCF 3 GC_bins 4 VCF_bins 5 VCF_counts
#use strict;
use warnings;

my %lengthHash;
my %GC_Hash;

open FAI, $ARGV[0] or die $!;

while(<FAI>)
{
	chomp;
	my @sps = split(/\t/);
	
	if(scalar @sps > 1)
	{
		$lengthHash{$sps[0]} = $sps[1];
	}
}

close FAI;

open GENOME, $ARGV[1] or die $!;

my $oldSC = "";
my $seq = "";

while(<GENOME>)
{
	chomp;
	my $subs = substr $_, 0, 1;
	
	if($subs eq ">")
	{
		if(length($seq) > 1)
		{
			my $len = length($seq);
			my $ads = 50;
			for(my $i = 0; $i < $len; $i += 50)
			{
				if($i > $len - 50)
				{
					$ads = $len - $i;
				}
				my $f50bp = substr $seq, $i, $ads;
				my $GC = ($f50bp =~ tr/GC//);
				$GC = $GC/$ads;
				$GC_Hash{$oldSC}{$i} = $GC;
			}
			$seq = "";
			print "Did GC $oldSC\n";
		}
		$oldSC = $_;
		$oldSC =~ tr/>//d;
	} else 
	{
		$seq .= $_;
	}
	
}

if(length($seq) > 1)
		{
			my $len = length($seq);
			my $ads = 50;
			for(my $i = 0; $i < $len; $i+=50)
			{
				if($i > $len - 50)
				{
					$ads = $len - $i;
				}
				my $f50bp = substr $seq, $i, $ads;
				my $GC = ($f50bp =~ tr/GC//);
				$GC = $GC/$ads;
				$GC_Hash{$oldSC}{$i} = $GC;
			}
			$seq = "";
			
			print "GC Did $oldSC\n";
		}

close GENOME;
my %wigHash;
my %countHash;

my $fiftyPlus = 0;
my %fiftyHash;
my $oldID  = "";

open VCF, $ARGV[2] or die $!;

while(<VCF>)
{
	chomp;
	
	if($. > 54)
	{
		my @sps = split(/\t/);
		my $len = length($sps[3]);
		my $len2 = length($sps[4]);
		my @sps2 = split(/_/, $sps[0]);
		
		if($len == $len2)
		{
			my $ct = 0; my $gt = 0; my $at = 0; my $cg = 0;
			
			for(my $i = 0; $i < $len; $i++)
			{
				my $test = (substr $sps[3], $i, 1) . (substr $sps[4], $i, 1);
				
				if($test eq "CT" || $test eq "GA" || $test eq "AG" || $test eq "TC")
				{
					$ct++;
				}
				elsif($test eq "GT" || $test eq "TG" || $test eq "CA" || $test eq "AC")
				{
					$gt++;
				}
				elsif($test eq "AT" || $test eq "TA")
				{
					$at++;
				}
				elsif($test eq "CG" || $test eq "GC")
				{
					$cg++;
				}
			}
			
			if($oldID eq $sps[0])
			{
				while($sps[1] > $fiftyPlus)
				{
					$fiftyPlus += 50;
					$fiftyHash{$sps2[0]}{$fiftyPlus}{"CT"} = 0;
					$fiftyHash{$sps2[0]}{$fiftyPlus}{"GT"} = 0;
					$fiftyHash{$sps2[0]}{$fiftyPlus}{"AT"} = 0;
					$fiftyHash{$sps2[0]}{$fiftyPlus}{"CG"} = 0;
				}
				
				$fiftyHash{$sps2[0]}{$fiftyPlus}{"CT"} += $ct;
				$fiftyHash{$sps2[0]}{$fiftyPlus}{"GT"} += $gt;
				$fiftyHash{$sps2[0]}{$fiftyPlus}{"AT"} += $at;
				$fiftyHash{$sps2[0]}{$fiftyPlus}{"CG"} += $cg;
				
			} else {
				print "Poly Doing " . $sps2[0] . "\n";
				$fiftyPlus = 50;
				$fiftyHash{$sps2[0]}{$fiftyPlus}{"CT"} = 0;
				$fiftyHash{$sps2[0]}{$fiftyPlus}{"GT"} = 0;
				$fiftyHash{$sps2[0]}{$fiftyPlus}{"AT"} = 0;
				$fiftyHash{$sps2[0]}{$fiftyPlus}{"CG"} = 0;
				
				print "Did " . $sps2[0] . "\n";
				while($sps[1] > $fiftyPlus)
				{
					$fiftyPlus += 50;
					$fiftyHash{$sps2[0]}{$fiftyPlus}{"CT"} = 0;
					$fiftyHash{$sps2[0]}{$fiftyPlus}{"GT"} = 0;
					$fiftyHash{$sps2[0]}{$fiftyPlus}{"AT"} = 0;
					$fiftyHash{$sps2[0]}{$fiftyPlus}{"CG"} = 0;
				}
				
				$fiftyHash{$sps2[0]}{$fiftyPlus}{"CT"} += $ct;
				$fiftyHash{$sps2[0]}{$fiftyPlus}{"GT"} += $gt;
				$fiftyHash{$sps2[0]}{$fiftyPlus}{"AT"} += $at;
				$fiftyHash{$sps2[0]}{$fiftyPlus}{"CG"} += $cg;
			}
			
			if(exists $countHash{$sps[0]})
			{
				$countHash{$sps[0]}{"CT"} += $ct;
				$countHash{$sps[0]}{"GT"} += $gt;
				$countHash{$sps[0]}{"AT"} += $at;
				$countHash{$sps[0]}{"CG"} += $cg;
			} else {
				$countHash{$sps[0]}{"CT"} = $ct;
				$countHash{$sps[0]}{"GT"} = $gt;
				$countHash{$sps[0]}{"AT"} = $at;
				$countHash{$sps[0]}{"CG"} = $cg;
			}
		}
		$oldID = $sps[0];
	}
	
}


close VCF;

my @keys = keys %GC_Hash;

open my $fout, '>', $ARGV[3] or die $!;

foreach(@keys)
{
	my $scaf = $_;
	
	my $fif = 0;
	
	while(exists $GC_Hash{$scaf}{$fif})
	{	
		print {$fout} $scaf . "\t" . $GC_Hash{$scaf}{$fif} . "\n";
		$fif += 50;
	}
}

close $fout;

open my $fout2, '>', $ARGV[4] or die $!;
my @keys2 = keys %fiftyHash;
print {$fout2} "Scaffold\tCT\tGT\tAT\tCG\n";

foreach(@keys2)
{
	my $scaf = $_;
	
	my $fif = 50;
	
	while(exists $fiftyHash{$scaf}{$fif}{"CG"})
	{	
		print {$fout2} $scaf . "\t" . $fiftyHash{$scaf}{$fif}{"CT"} . "\t" . $fiftyHash{$scaf}{$fif}{"GT"} . "\t" . $fiftyHash{$scaf}{$fif}{"AT"} . "\t" . $fiftyHash{$scaf}{$fif}{"CG"} . "\n";
		$fif += 50;
	}
}

close $fout2;

open my $countOut, '>', $ARGV[5] or die $!;
my @keys3 = keys %countHash;
print {$countOut} "Scaffold\tLength\tCT\tGT\tAT\tCG\n";

foreach(@keys3)
{
	my $scaf = $_;
	my $cnt = $countHash{$scaf}{"CT"} + $countHash{$scaf}{"GT"} + $countHash{$scaf}{"AT"} + $countHash{$scaf}{"CG"};
	
	my $len = $lengthHash{$scaf};
	my $cto = ($countHash{$scaf}{"CT"} / $cnt) * 100;
	my $gto = ($countHash{$scaf}{"GT"} / $cnt) * 100;
	my $ato = ($countHash{$scaf}{"AT"} / $cnt) * 100;
	my $cgo = ($countHash{$scaf}{"CG"} / $cnt) * 100;
	
	print {$countOut} $scaf . "\t" . $len  . "\t" . $cto . "\t" . $gto . "\t" . $ato . "\t" . $cgo ."\n";
}

close $countOut;
