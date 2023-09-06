#!/usr/bin/perl -w
# made by Konstantin V. Gunbin aka GenKVG
# required OS: Linux
# required tools (installed in /usr/local/bin by default):
#	1) pigz - https://zlib.net/pigz/
# required script location: the directory containing VCF files
# run: ./step0.pl

while ($file=<*.vcf.gz>)
{
print "$file\n";
$mnp=$file;
$mnp=~s/\.gz$/\.mnp/;
open (O, ">$mnp");
open (V, "zcat $file |");
while (<V>)
{
if ($_=~/^#/) {print O "$_"}
else 
{
chomp;
@V=split(/\t/);
$chr=$V[0];
@V3=split(//, $V[3]);
@V4=split(//, $V[4]);
$tail=""; for ($i=5;$i<$#V+1;$i++) {$tail.="\t$V[$i]"}
if ($V[6] eq "PASS")
{
if ($#V3>0 and $#V3==$#V4)
{
for ($i=0;$i<$#V3+1;$i++)
    {
    $posn=$V[1]+$i;
    print O "$chr\t$posn\t\.\t$V3[$i]\t$V4[$i]$tail\n";
    }
}
else {print O "$_\n"}
}
}
}
close O;
system ("pigz $mnp");
}
