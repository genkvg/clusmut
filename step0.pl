#!/usr/bin/perl -w
# made by Konstantin V. Gunbin aka GenKVG
# required OS: Linux
# required tools (installed in /usr/local/bin by default):
#	1) convert2bed - https://bedops.readthedocs.io/en/latest/index.html
#	2) sambamba - https://lomereiter.github.io/sambamba/
#	3) jvarkit [Biostar322664] - http://lindenb.github.io/jvarkit/JvarkitCentral.html [http://lindenb.github.io/jvarkit/Biostar322664.html]  
# required script location: the directory containing VCF and BAM files
# run: ./step0.pl

%vcf=();
open (L, "ls -1 *.vcf |");
while (<L>)
{
chomp;
@L=split(/!SYMB!/);		#set the combination of symbols (!SYMB!) to parse VCF file name into the array, the first element of array must be tumor ID 
$vcf{$L[0]}=$_;
}
close L;

foreach $v (sort {$a cmp $b} keys %vcf) 
{ 
$bed=$vcf{$v};
$bed=~s/\.vcf$/\.bed/; 
$bam=$v."!BFPF!";		#BAM files postfix (!BFPF!)
unless (-e "$bam.Smutper_read_nothr_srt_supp.bam")
{
unless (-e "$bed") {system ("convert2bed --input=vcf --output=bed < $vcf{$v} > $bed")}
system ("sambamba view -t 20 -F \"([NM] >= 1)\" -f bam -L $bed -o $bam.Smutper_read_nothr.bam $bam.bam"); 
system ("mkdir tmp");
system ("sambamba sort -m 8G -t 20 --tmpdir=tmp -o $bam.Smutper_read_nothr_srt.bam $bam.Smutper_read_nothr.bam"); 
system ("java -jar /usr/local/bin/biostar322664.jar -nm -V $vcf{$v} --samoutputformat bam -o $bam.Smutper_read_nothr_srt_supp.bam $bam.Smutper_read_nothr_srt.bam");
system ("rm $bam.Smutper_read_nothr_srt.bam* $bam.Smutper_read_nothr.bam*");
system ("rm -rf tmp");
}
}
