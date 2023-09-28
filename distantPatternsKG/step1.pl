#!/usr/bin/perl -w
# made by Konstantin V. Gunbin aka GenKVG
# required OS: Linux
# required tools (installed in /usr/local/bin by default):
#	1) pigz - https://zlib.net/pigz/
# required script location: the directory containing tabix indexes [see step0.pl]
# run: ./step1.pl

while ($file=<*.gz>)
{
$chr=$file;
$chr=~s/\.gz//;

%trinuc=();
$n=0;
open (I, "zcat $file |");
while (<I>)
{
chomp;
$n++;
print STDERR "$chr - $n\r";
@I=split(/\t/);
#8       10002   10002   G       C       A
$trinuc{"$I[3]$I[4]$I[5]"}{$I[1]}=1;
}
close I;

system ("mkdir trinucs");

foreach $tn (keys %trinuc)
{
open (O, "+>>trinucs/$tn.lst");
foreach $pos (sort {$a <=> $b} keys %{$trinuc{$tn}})
{
print O "$chr\t$pos\n";
}
close O;
}
}

system ("pigz -p8 *.lst");
