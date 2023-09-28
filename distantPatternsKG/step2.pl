#!/usr/bin/perl -w
# made by Konstantin V. Gunbin aka GenKVG
# required OS: Linux
# required RAM: 200Gb of free RAM
# required tools (installed in /usr/local/bin by default):
#	1) Math::Random::MT::Auto Perl module - https://metacpan.org/pod/Math::Random::MT::Auto
#	2) Mappable position BED3 file - e.g. https://genome.ucsc.edu/cgi-bin/hgTables?db=mm9&hgta_group=map&hgta_track=wgEncodeMapability&hgta_table=wgEncodeCrgMapabilityAlign75mer
#	3) toPERM.tsv file [see step2.pl in observedVSexpected directory]
#	4) pigz - https://zlib.net/pigz/
# required script location: the directory containing '*.lst.gz' files [see step1.pl]
# run: ./step2.pl Mappable_position.BED3 toPERM.tsv distance_from_real_mutations[e.g. 10000]

use Math::Random::MT::Auto qw(shuffle srand :!auto);

print STDERR "Real Positions...\n";
%posrealtrinuc=();
open (I, "$ARGV[1]");
while (<I>)
{
chomp;
@I=split(/\t/);
$hundK=int($I[1]/500000);
$posrealtrinuc{$I[3]}{$I[0]}{$hundK}{$I[1]}=1;
}
close I;
print STDERR "Real Positions Hashed!\n";

print STDERR "Mappable Positions...\n";
%noblacklist=();
open (I, "$ARGV[0]");
while (<I>)
{
chomp;
@I=split(/\t/);
print STDERR "$I[0]\r";
for ($i=$I[1];$i<$I[2]+1;$i++)
    {
    $noblacklist{$I[0]}{$i}=1;
    }
}
close I;
print STDERR "Mappable Positions Hashed!\n";

while ($file=<*.lst.gz>)
{
@mutations=();
@trinuc=split(/\./, $file);
$vcf=$file;
$vcf=~s/\.lst\.gz/\.vcf/;
print STDERR "$file\n";
open (I, "zcat $file |");
@patt=split(/\./, $file);
@P=split(//, $patt[0]);
while (<I>)
{
chomp;
@II=split(/\t/);
$hundK=int($II[1]/500000);
if (exists $noblacklist{$II[0]}{$II[1]} and exists $posrealtrinuc{$trinuc[0]}{$II[0]}{$hundK})
    {
    $min=500000;
    foreach $ps (keys %{$posrealtrinuc{$trinuc[0]}{$II[0]}{$hundK}})
	{
	$nim=abs($ps-$II[1]);
	if ($min>$nim) {$min=$nim}
	}
    if ($min<=$ARGV[2]) {push @mutations, "$_";}
    }
}
close I;

@N=("A","T","G","C");
@M=(); foreach $n (@N) {unless ($n eq $P[1]) {push @M, $n}}

for ($s=0;$s<10;$s++) {shuffle(\@mutations)}

open (O, ">$vcf");
print O qq(##fileformat=VCFv4.1
##reference=Homo_sapiens.GRCh37.GATK.illumina.fasta
##source=$file
##fileDate=29.07.2022
##FILTER=<ID=Default,Description="Default">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##contig=<ID=1,length=249250621,assembly="hg19">
##contig=<ID=2,length=243199373,assembly="hg19">
##contig=<ID=3,length=198022430,assembly="hg19">
##contig=<ID=4,length=191154276,assembly="hg19">
##contig=<ID=5,length=180915260,assembly="hg19">
##contig=<ID=6,length=171115067,assembly="hg19">
##contig=<ID=7,length=159138663,assembly="hg19">
##contig=<ID=8,length=146364022,assembly="hg19">
##contig=<ID=9,length=141213431,assembly="hg19">
##contig=<ID=10,length=135534747,assembly="hg19">
##contig=<ID=11,length=135006516,assembly="hg19">
##contig=<ID=12,length=133851895,assembly="hg19">
##contig=<ID=13,length=115169878,assembly="hg19">
##contig=<ID=14,length=107349540,assembly="hg19">
##contig=<ID=15,length=102531392,assembly="hg19">
##contig=<ID=16,length=90354753,assembly="hg19">
##contig=<ID=17,length=81195210,assembly="hg19">
##contig=<ID=18,length=78077248,assembly="hg19">
##contig=<ID=19,length=59128983,assembly="hg19">
##contig=<ID=20,length=63025520,assembly="hg19">
##contig=<ID=21,length=48129895,assembly="hg19">
##contig=<ID=22,length=51304566,assembly="hg19">
##contig=<ID=X,length=155270560,assembly="hg19">
##contig=<ID=Y,length=59373566,assembly="hg19">
##contig=<ID=MT,length=16569,assembly="hg19">
);
$n=-1;
print O "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSMP\n";
%mutations1M=();
for ($i=0;$i<300000;$i++) {@MT=split(/\t/,$mutations[$i]); $mutations1M{$MT[0]}{$MT[1]}=1}
foreach $m (sort {$a cmp $b} keys %mutations1M)
    {
    foreach $s (sort {$a <=> $b} keys %{$mutations1M{$m}})
    {
    $n++;
    if ($n>2) {$n=-1; $n++}
    print O "$m\t$s\t.\t$P[1]\t$M[$n]\t100\tPASS\tAC=1;AF=0.5;DP=100\tGT:AD:DP\t0/1:50,50:100\n";
    }
    }
close O;
}

system ("pigz -p8 *.vcf");
