#!/usr/bin/perl -w
# made by Konstantin V. Gunbin aka GenKVG
# required OS: Linux
# required tools (installed in /usr/local/bin by default):
#	1) tabix - http://www.htslib.org/download/
#	2) samtools - http://www.htslib.org/download/
#	3) human reference genome [hg19 or compatible version] - https://gatk.broadinstitute.org/hc/en-us/articles/360035890711-GRCh37-hg19-b37-humanG1Kv37-Human-Reference-Discrepancies
# required script location: empty directory
# run: ./step0.pl

%chrs=();
open (I, "/path/to/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai");			#the directory with human reference genome (samtools *fai index) file must be set up!
while (<I>)
{
chomp;
@I=split(/\t/);
$chrs{$I[0]}=$I[1];
}
close I;
print STDERR "chrs length\n";

@chrs=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y);
#@chrs=(MT);

foreach $chr (@chrs)
{
print STDERR "#########################$chr#########################\n";
$CHR=" ";
$l=0;
open (C, "samtools faidx /path/to/Homo_sapiens.GRCh37.GATK.illumina.fasta $chr |");	#the directory with human reference genome file must be set up!
while (<C>)
{
$l++;
print STDERR "$l\r";
chomp;
unless ($_=~/^>/) {$CHR.=$_}
}
close C;

open O, '|-', "bgzip -c > $chr.gz";
for ($i=2;$i<$chrs{$chr}+1;$i++)
{
$pattern=substr $CHR, $i-1, 3;
unless ($pattern=~/N/)
{
@pattern=split(//, $pattern);
print O "$chr\t$i\t$i\t$pattern[0]\t$pattern[1]\t$pattern[2]\n";
}
}
close O;
system ("tabix -s 1 -b 2 -e 3 $chr.gz");
}
