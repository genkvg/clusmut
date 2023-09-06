#!/home/genkvg/perl/bin/perl -w
# made by Konstantin V. Gunbin aka GenKVG
# required OS: Linux
# required tools (installed in /usr/local/bin by default):
#	1) bedtools - https://bedtools.readthedocs.io/en/latest/
#	2) human reference genome [hg19 or compatible version] - https://gatk.broadinstitute.org/hc/en-us/articles/360035890711-GRCh37-hg19-b37-humanG1Kv37-Human-Reference-Discrepancies
# required script location: the directory containing filtered VCF files [see step0.pl]
# run: ./step1.pl

open (O, ">toBED.bed");
while ($file=<*.vcf.mnp.gz>)
{
@sample=split(/!SYMB!/,$file);		#set the combination of symbols (!SYMB!) to parse BED file name into the array, the first element of array must be tumor ID
print STDERR "$file\n";
open (V, "zcat $file |");
while (<V>)
{
unless ($_=~/^#/)
{
chomp;
@V=split(/\t/);
$chr=$V[0];
$pos=$V[1];
$pre=$pos-2;
$post=$pos+1;
$ref=$V[3];
@ref=split(//,$ref);
$alt=$V[4];
@alt=split(//,$alt);

if ($#ref==0 and $#ref==$#alt)
{
print O "$chr\t$pre\t$post\t$sample[0]\t$ref\t$alt\n";
}
}
}
}
close O;

sleep(15);
system ("bedtools getfasta -fi Homo_sapiens.GRCh37.GATK.illumina.fasta -bed toBED.bed > toBED.fa");
