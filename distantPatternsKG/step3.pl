#!/home/genkvg/perl/bin/perl -w
# made by Konstantin V. Gunbin aka GenKVG
# required OS: Linux
# required tools (installed in /usr/local/bin by default):
#	1) human reference genome [hg19 or compatible version] - https://gatk.broadinstitute.org/hc/en-us/articles/360035890711-GRCh37-hg19-b37-humanG1Kv37-Human-Reference-Discrepancies
# required script location: the directory containing filtered 'vcf.gz' files [see step2.pl]
# run: ./step3.pl max_distance_to_dinucleotide[e.g. 500] OUTDISTANCES.random

@dinuc=("CC","TT","AA","GG","CT","TC","GA","AG");

open (O, ">$ARGV[1].$ARGV[0].FCON");
print O "Dinucleotide\tChange";
for ($i=-$ARGV[0]+1;$i<$ARGV[0];$i++) {print O "\t$i"}
print O "\n";

%patts=();
%distTOpirdim=();
%changes=();
while ($file=<*.vcf.gz>)
{
@F=();
if ($file=~/\.vcf/) {@F=split(/\.vcf/,$file)}
$smp=$F[0];

print STDERR "$smp\n";
open (B0, ">$smp.0.bed");
open (V, "zcat $file |");
while (<V>)
{
unless ($_=~/^#/)
{
chomp;
@V=split(/\t/);
$chr=$V[0];
$pos=$V[1];
@V3=split(//, $V[3]);
@V4=split(//, $V[4]);
if ($#V3==0 and $#V3==$#V4)
{unless ($chr eq "MT"){
$pre=$pos-2;
$post=$pos+1;
print B0 "$chr\t$pre\t$post\n";
}}
}
}
close V;
close B0;

system ("bedtools getfasta -fi Homo_sapiens.GRCh37.GATK.illumina.fasta -bed $smp.0.bed > $smp.0.fa");

open (S, "$smp.0.fa");
while (<S>)
{
if ($_=~/^>([0-9XY]+):(\d+)-/) {$ch=$1;$ps=$2;$ps=$ps+2}
else
    {
    chomp;
    $patts{"$ch:$ps"}=$_;
    }
}
close S;

open (B1, ">$smp.1.bed");
open (B2, ">$smp.2.bed");
open (V, "zcat $file |");
while (<V>)
{
unless ($_=~/^#/)
{
chomp;
@V=split(/\t/);
$chr=$V[0];
unless ($chr eq "MT")
{
$pos=$V[1];
@V3=split(//, $V[3]);
@V4=split(//, $V[4]);
if ($#V3==0 and $#V3==$#V4)
{
unless (exists $patts{"$chr:$pos"}){print STDERR "$chr:$pos\n"}
@P=split(//, $patts{"$chr:$pos"});
if ($V[3] eq $P[1])
{
$change="$P[0]\($P[1]>$V[4]\)$P[2]";
$posF=$pos;
$end=$posF+$ARGV[0];

print B1 "$chr\t$posF\t$end\t$change\n";

$posB=$pos-1;
$end=$posB-$ARGV[0];

unless ($end<0) {print B2 "$chr\t$end\t$posB\t$change\n"}
}
}
}
}
}
close V;
close B1;
close B2;

system ("bedtools getfasta -fi Homo_sapiens.GRCh37.GATK.illumina.fasta -name -bed $smp.1.bed > $smp.1.fa");
system ("bedtools getfasta -fi Homo_sapiens.GRCh37.GATK.illumina.fasta -name -bed $smp.2.bed > $smp.2.fa");


open (S, "$smp.1.fa");
while (<S>)
{
if ($_=~/^>([ATGC]\([ATGC]>[ATGC]\)[ATGC])::/) {$change="$1"; $changes{$change}=1; $changes{all}=1}
else
{
chomp;
$seq=$_;
foreach $di (@dinuc)
{
$pk = index $seq, $di;
$distTOpirdim{$di}{$pk}{$change}+=1;
$distTOpirdim{$di}{$pk}{all}+=1;
$distTOpirdim{all}+=1;
}
}
}
close S;

open (S, "$smp.2.fa");
while (<S>)
{
if ($_=~/^>([ATGC]\([ATGC]>[ATGC]\)[ATGC])::/) {$change="$1"; $changes{$change}=1; $changes{all}=1}
else 
{
chomp;
$seq=$_;
@seq=split(//);
@revseq=reverse(@seq);
$revseq=join("", @revseq);
foreach $di (@dinuc)
{
$pk = index $revseq, $di;
$distTOpirdim{$di}{"-$pk"}{$change}+=1;
$distTOpirdim{$di}{"-$pk"}{all}+=1;
$distTOpirdim{all}+=1;
}
}
}
close S;

system ("rm $smp.0.bed $smp.1.bed $smp.2.bed $smp.1.fa $smp.2.fa $smp.0.fa");
}

foreach $di (@dinuc)
{
foreach $ch (sort {$a cmp $b} keys %changes)
{
print O "$di\t$ch";
for ($i=-$ARGV[0]+1;$i<$ARGV[0];$i++)
{
if (exists $distTOpirdim{$di}{$i}{$ch}) {$pc=sprintf("%.0f", $distTOpirdim{$di}{$i}{$ch}); print O "\t$pc"} else {print O "\t0.000"}
}
print O "\n";
}
}
close O;
