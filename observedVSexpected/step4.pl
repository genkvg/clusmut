#!/home/genkvg/perl/bin/perl -w
# made by Konstantin V. Gunbin aka GenKVG
# required OS: Linux
# required tools (installed in /usr/local/bin by default):
#	1) Math::Random::MT::Auto Perl module - https://metacpan.org/pod/Math::Random::MT::Auto
#	2) Parallel::ForkManager Perl module - https://metacpan.org/pod/Parallel::ForkManager  
# required script location: the directory containing toPERM.tsv [see step2.pl] and 'patterns' files
# run: ./step4.pl toPERM.tsv

use Math::Random::MT::Auto qw(shuffle srand :!auto);
use Parallel::ForkManager;

%patterns=();
open (I, "patterns");
while (<I>)
{
chomp;
@I=split(/\t/);
$patterns{"$I[0]\t$I[1]\t$I[2]\t$I[3]"}=$I[4];
}
close I;

@mutations=();
@samples=();
open (I, "$ARGV[0]");
while (<I>)
{
chomp;
@I=split(/\t/);
push @mutations, "$I[0]\t$I[1]\t$I[3]---$I[4]";
push @samples, "$I[2]";
}
close I;

my $pm = new Parallel::ForkManager(10);
for ($it=0;$it<$ARGV[1];$it++)
    {
    $pm->start and next;
    srand();
    &run();
    $pm->finish;
    }
$pm->wait_all_children;

sub run
{
my(@ssamples, $s, %mutationsSamples, $q, @mut, %patternMatrix, $samp, $ch, @samplechr, $pos, $i, @A, @B, $j, $dist, $key, $o);
@ssamples=@samples;
for ($s=0;$s<10;$s++) {shuffle(\@ssamples)}

%mutationsSamples=();
for ($q=0;$q<$#mutations+1;$q++)
    {
    @mut=split(/\t/, $mutations[$q]);
    $mutationsSamples{$ssamples[$q]}{$mut[0]}{$mut[1]}=$mut[2];
    }
%patternMatrix=();
foreach $samp (keys %mutationsSamples)
    {
    foreach $ch (sort {$a cmp $b} keys %{$mutationsSamples{$samp}})
	{
	@samplechr=();
	foreach $pos (sort {$a <=> $b} keys %{$mutationsSamples{$samp}{$ch}})
	    {
	    push @samplechr, "$pos\t$mutationsSamples{$samp}{$ch}{$pos}";
	    }
	for ($i=0;$i<$#samplechr+1;$i++)
	    {
	    @A=split(/\t/, $samplechr[$i]);
	    for ($j=$i+1;$j<$#samplechr+1;$j++)
		{
		@B=split(/\t/, $samplechr[$j]);
		$dist=$B[0]-$A[0];
		if ($dist>100) {last}
		else
		    {
		    $A[1]=~s/---/\t/;
		    $B[1]=~s/---/\t/;
		    $patternMatrix{"$A[1]\t$B[1]"}{$dist}+=1;
		    }
		}
	    }
	}
    }

open (O, ">rnd.$it");
foreach $key (sort {$patterns{$a} <=> $patterns{$b}} keys %patterns)
    {
    print O "$key\t$patterns{$key}";
    for ($o=0;$o<101;$o++)
	{
	if (exists $patternMatrix{$key}{$o}) {print O "\t$patternMatrix{$key}{$o}"} else {print O "\t0"}
	}
    print O "\n";
    }
close O;
}
