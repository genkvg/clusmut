#!/usr/bin/perl -w
# made by Konstantin V. Gunbin aka GenKVG
# required OS: Linux
# required tools (installed in /usr/local/bin by default):
#	1) samtools - http://www.htslib.org/download/
#	2) human reference genome [hg19 or compatible version] - https://gatk.broadinstitute.org/hc/en-us/articles/360035890711-GRCh37-hg19-b37-humanG1Kv37-Human-Reference-Discrepancies
# required script location: the directory containing BED and BAM (original and *supp [see step0.pl]) files
# run: ./step1.pl

%bed=();
open (L, "ls -1 *.bed |");
while (<L>)
{
chomp;
@L=split(/!SYMB!/);				#set the combination of symbols (!SYMB!) to parse BED file name into the array, the first element of array must be tumor ID
$bed{$L[0]}=$_;
}
close L;

%bamsOri=();
open (L, "ls -1 *!OBFP!.bam |");		#original BAM files postfix (!OBFP!)
while (<L>)
{
chomp;
@L=split(/!SYMB!/);				#set the combination of symbols (!SYMB!) to parse original BAM file name into the array, the first element of array must be tumor ID
$bamsOri{$L[0]}=$_;
}
close L;

%bamsAlt=();
open (L, "ls -1 *supp.bam |");
while (<L>)
{
chomp;
@L=split(/!SYMB!/);				#set the combination of symbols (!SYMB!) to parse *supp BAM file name into the array, the first element of array must be tumor ID
push @{$bamsAlt{$L[0]}}, $_;
}
close L;

%table=();

foreach $v (sort {$a cmp $b} keys %bed) 
{ 
$TID=$v;
if (exists $bamsOri{$v} and $bamsAlt{$v})
{
$oribam=$bamsOri{$v};
foreach $altbam (@{$bamsAlt{$v}})
{
print STDERR "$altbam\n";
print "$oribam\n";
$altmpileup=$distances=$altbam;
$altmpileup=~s/\.bam$/\.mpileup/;
$distances=~s/\.bam$/\.bimut.dists/;
$orimpileup=$oribam;
$orimpileup=~s/\.bam$/\.mpileup/;

unless (-e "$altmpileup")
{
system ("samtools mpileup -f Homo_sapiens.GRCh37.GATK.illumina.fasta -l $bed{$v} --output-QNAME $altbam -o $altmpileup");
}
$linepre=""; $chrpre=""; $pospre=""; @allelespre=(); @readspre=(); $refpre="";
%dists=();
$filetype="alt";
$mpileup=$altmpileup;
&readsCount();


open (O, ">$distances");
print O "Distance\tBoth reads #\t5' reads #\t3' reads #\n";
foreach $dist (sort {$a<=>$b} keys %dists)
{
print O "$dist\t$dists{$dist}{both}\t$dists{$dist}{5}\t$dists{$dist}{3}\n";
}
close O;
}

unless (-e "$orimpileup")
{
system ("samtools mpileup -f Homo_sapiens.GRCh37.GATK.illumina.fasta -l $bed{$v} --output-QNAME $oribam -o $orimpileup");
}
$linepre=""; $chrpre=""; $pospre=""; @allelespre=(); @readspre=(); $refpre="";
$filetype="ori";
$mpileup=$orimpileup;
&readsCount();
}
}

open (O, ">table");
print O "TID\tchr\tpos 5'\tref 5'\talt 5'\tpos 3'\tref 3'\talt 3'\tdistance\tAll reads #\tBimutated reads #\t5' mutated reads #\t3' mutated reads #\n";
foreach $tid (sort {$a cmp $b} keys %table)
{
foreach $chr (sort {$a cmp $b} keys %{$table{$tid}})
{
foreach $poss (sort {$a cmp $b} keys %{$table{$tid}{$chr}})
{
@poss=split(/\t/,$poss);
print O "$tid\t$chr\t$poss[0]\t$table{$tid}{$chr}{$poss}{r5}\t$table{$tid}{$chr}{$poss}{a5}\t$poss[1]\t$table{$tid}{$chr}{$poss}{r3}\t$table{$tid}{$chr}{$poss}{a3}\t$table{$tid}{$chr}{$poss}{d}\t$table{$tid}{$chr}{$poss}{a}\t$table{$tid}{$chr}{$poss}{b}\t$table{$tid}{$chr}{$poss}{5}\t$table{$tid}{$chr}{$poss}{3}\n";
}
}
}
close O;


sub readsCount
{
open (M, "$mpileup");
while (<M>)
{
#1       1059398 T       2       AA      DB      V300090768L2C004R0520447184,V300090768L3C002R0020194945
chomp;
@M=split(/\t/);
if ($chrpre eq "") 
    {
    $linepre=$_;
    $chrpre=$M[0]; $pospre=$M[1]; $M[2]=~tr/a-z/A-Z/; $refpre=$M[2];
    $M[4]=~tr/a-z/A-Z/; $M[4]=~s/\^.{1}|\$//g; 
    while ($M[4]=~/([\+\-]\d+)/) {$d=$e=$1; $d="\\$d"; $e=~s/[\+\-]//; print STDERR "$M[4] ||"; $M[4]=~s/$d[A-Z]{$e}//; print STDERR "$M[4] $d\n"}
    @allelespre=split(//, $M[4]);
    @readspre=split(/,/, $M[6]);
    }
else
    {
    $line=$_;
    $chr=$M[0]; $pos=$M[1]; $M[2]=~tr/a-z/A-Z/; $ref=$M[2];
    $M[4]=~tr/a-z/A-Z/; $M[4]=~s/\^.{1}|\$//g; 
    while ($M[4]=~/([\+\-]\d+)/) {$d=$e=$1; $d="\\$d"; $e=~s/[\+\-]//; print STDERR "$M[4] ||"; $M[4]=~s/$d[A-Z]{$e}//; print STDERR "$M[4] $d\n"}
    @alleles=split(//, $M[4]);
    @reads=split(/,/, $M[6]);
    if ($chrpre eq $chr and $pos-$pospre<101)
	{ if ($#alleles==$#reads and $#allelespre==$#readspre) {
	$dist=$pos-$pospre;
	%end5reads=();
	%end5readsref=();
	%allelespre=();
	for ($i=0;$i<$#allelespre+1;$i++)
	    {
	    if ($allelespre[$i]=~/[A-Z]/) {if ($allelespre[$i] ne $refpre) {
		$allelespre{$allelespre[$i]}+=1;
		$end5reads{$readspre[$i]}=1;
		}}
	    else {$end5readsref{$readspre[$i]}=1}
	    }
	%end3reads=();
	%end3readsref=();
	%alleles=();
	for ($i=0;$i<$#alleles+1;$i++)
	    {
	    if ($alleles[$i]=~/[A-Z]/) {if ($alleles[$i] ne $ref) {
		$alleles{$alleles[$i]}+=1;
		$end3reads{$reads[$i]}=1;
		}}
	    else {$end3readsref{$reads[$i]}=1}
	    }
	@allelespre=sort{$allelespre{$b}<=>$allelespre{$a}} keys %allelespre; #$allelespre=join(".", @allelespre);
	@alleles=sort{$alleles{$b}<=>$alleles{$a}} keys %alleles; #$alleles=join(".", @alleles);
	$both=0;$end5=0;$end3=0;$bothref=0;
	foreach $r (keys %end5reads)
	    {if (exists $end3reads{$r}) {$both+=1}}
	foreach $r (keys %end5reads)
	    {if (exists $end3readsref{$r}) {$end5+=1}}
	foreach $r (keys %end5readsref)
	    {if (exists $end3reads{$r}) {$end3+=1}}
	foreach $r (keys %end5readsref)
	    {if (exists $end3readsref{$r}) {$bothref+=1}}
	if ($filetype eq "alt")
	{
	$dists{$dist}{both}+=$both;
	$dists{$dist}{5}+=$end5;
	$dists{$dist}{3}+=$end3;
	$table{$TID}{$chr}{"$pospre\t$pos"}{5}=$end5;
	$table{$TID}{$chr}{"$pospre\t$pos"}{"r5"}=$refpre;
	$table{$TID}{$chr}{"$pospre\t$pos"}{"a5"}=$allelespre[0];
	$table{$TID}{$chr}{"$pospre\t$pos"}{3}=$end3;
	$table{$TID}{$chr}{"$pospre\t$pos"}{"r3"}=$ref;
	$table{$TID}{$chr}{"$pospre\t$pos"}{"a3"}=$alleles[0];
	$table{$TID}{$chr}{"$pospre\t$pos"}{b}=$both;
	}
	elsif ($filetype eq "ori")
	{
	$table{$TID}{$chr}{"$pospre\t$pos"}{a}=$bothref+$both+$end5+$end3;
	$table{$TID}{$chr}{"$pospre\t$pos"}{d}=$dist;
	}
	}
	else {print STDERR "$#allelespre $#readspre -- $linepre\n$#alleles $#reads -- $line\n"}
	}
    $linepre=$_;
    $chrpre=$M[0]; $pospre=$M[1]; $M[2]=~tr/a-z/A-Z/; $refpre=$M[2];
    $M[4]=~tr/a-z/A-Z/; $M[4]=~s/\^.{1}|\$//g; 
    while ($M[4]=~/([\+\-]\d+)/) {$d=$e=$1; $d="\\$d"; $e=~s/[\+\-]//; print STDERR "$M[4] ||"; $M[4]=~s/$d[A-Z]{$e}//; print STDERR "$M[4] $d\n"}
    @allelespre=split(//, $M[4]);
    @readspre=split(/,/, $M[6]);
    }
}
close M;
}