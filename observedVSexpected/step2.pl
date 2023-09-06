#!/usr/bin/perl -w
# made by Konstantin V. Gunbin aka GenKVG
# required OS: Linux
# required tools (installed in /usr/local/bin by default): no requirenments
# required script location: the directory containing toBED.bed and toBED.fa files [see step1.pl]
# run: ./step2.pl toBED.fa toBED.bed 1> toPERM.tsv

%trios=();
$l=0;
open (I, "$ARGV[0]");
while (<I>)
{
chomp;
if ($_=~/^>([0-9XY]+):(\d+)-(\d+)$/)
    {
    $c=$1; $s=$2; $e=$3;
    $key="$c\t$s\t$e";
    }
else
    {
    $l++;
    print STDERR "$l\r";
    @I=split(//);
    @{$trios{$key}}=@I;
    }
}
close I;
print STDERR "\n";

$l=$k=0;
open (I, "$ARGV[1]");
while (<I>)
{
$l++;
print STDERR "$l";
chomp;
@I=split(/\t/);
$key="$I[0]\t$I[1]\t$I[2]";
if (exists $trios{$key})
    {
    if ($trios{$key}[1] eq $I[4])
	{
	$k++;
	print STDERR " - $k\r";
	$ref=join("", @{$trios{$key}});
	$alt="$trios{$key}[0]$I[5]$trios{$key}[2]";
	$s=$I[1]+2;
	print "$I[0]\t$s\t$I[3]\t$ref\t$alt\n";
	}
    }
}
close I;
print STDERR "\n";

