#!/usr/bin/perl -w
# made by Konstantin V. Gunbin aka GenKVG
# required OS: Linux
# required tools (installed in /usr/local/bin by default):
#	1) tabix - http://www.htslib.org/download/
#	2) Bio::DB::HTS::Tabix Perl module - https://metacpan.org/pod/Bio::DB::HTS::Tabix
#	3) tabix index of human reference genome containing [see step2.pl]  
# required script location: the directory containing 'table' file
# run: ./step3.pl table 1> table.differentReadMut 2> table.sameReadMut

use Bio::DB::HTS::Tabix;

%bakhyTableMono=%bakhyTableBi=();
open (F, "$ARGV[0]");
while (<F>)
{
$l++;
if ($l>1)
{
chomp;
@F=split(/\t/);
if ($#F==12)
    {
    $z=$F[10]/($F[10]+$F[11]+$F[12]);
    if ($F[1]=~/[0-9a-zA-Z]/ and $F[4]=~/[ATGC]/ and $F[7]=~/[ATGC]/ and $F[3]=~/[ATGC]/ and $F[6]=~/[ATGC]/)
    {
    $queryA="$F[1]:$F[2]-$F[2]";
    $queryB="$F[1]:$F[5]-$F[5]";
    $tabix = Bio::DB::HTS::Tabix->new(filename => "/dir/with/posTabixGRCh37/$F[1].gz", use_tmp_dir => 1); #the directory with tabix indexed files containing trinucleotide positions must be set up!
    $iter = $tabix->query("$queryA");
    $pattA="";
    if (defined $iter)
        {
        @R=(); while(my $R = $iter->next) {push @R, $R}
        foreach $r (@R) {@r=split(/\t/, $r); if ($r[4] eq $F[3]) {$pattA="$r[3]$r[4]$r[5]"; $pattAm="$r[3]$F[4]$r[5]"}}
        }
    @pattA=split(//, $pattA);
    $iter = $tabix->query("$queryB");
    $pattB="";
    if (defined $iter)
        {
        @R=(); while(my $R = $iter->next) {push @R, $R}
        foreach $r (@R) {@r=split(/\t/, $r); if ($r[4] eq $F[6]) {$pattB="$r[3]$r[4]$r[5]"; $pattBm="$r[3]$F[7]$r[5]"}}
        }
    @pattB=split(//, $pattB);
    if ($#pattA==2 and $#pattB==2)
        {
#        print STDERR "$l\r";
        if ($z==0) {$bakhyTableMono{$pattA}{$pattAm}{$pattB}{$pattBm}{$F[8]}+=1}
        if ($z>0) {$bakhyTableBi{$pattA}{$pattAm}{$pattB}{$pattBm}{$F[8]}+=1}
	#print STDERR "$pattA>$pattAm $pattB>$pattBm $F[8]\n";
	}
    }
    }
}
}
close F;

$count=0;
foreach $A (sort {$a cmp $b} keys %bakhyTableMono)
{
foreach $Am (sort {$a cmp $b} keys %{$bakhyTableMono{$A}})
{
foreach $B (sort {$a cmp $b} keys %{$bakhyTableMono{$A}{$Am}})
{
foreach $Bm (sort {$a cmp $b} keys %{$bakhyTableMono{$A}{$Am}{$B}})
{
$count++;
print "$A\t$Am\t$B\t$Bm\t$count";
$string="";
for ($i=0;$i<101;$i++)
    {
    if (exists $bakhyTableMono{$A}{$Am}{$B}{$Bm}{$i}) {$string.="\t$bakhyTableMono{$A}{$Am}{$B}{$Bm}{$i}"}
    else {$string.="\t0"}
    }
print "$string\n"
}
}
}
}

$count=0;
foreach $A (sort {$a cmp $b} keys %bakhyTableBi)
{
foreach $Am (sort {$a cmp $b} keys %{$bakhyTableBi{$A}})
{
foreach $B (sort {$a cmp $b} keys %{$bakhyTableBi{$A}{$Am}})
{
foreach $Bm (sort {$a cmp $b} keys %{$bakhyTableBi{$A}{$Am}{$B}})
{
$count++;
print STDERR "$A\t$Am\t$B\t$Bm\t$count";
$string="";
for ($i=0;$i<101;$i++)
    {
    if (exists $bakhyTableBi{$A}{$Am}{$B}{$Bm}{$i}) {$string.="\t$bakhyTableBi{$A}{$Am}{$B}{$Bm}{$i}"}
    else {$string.="\t0"}
    }
print STDERR "$string\n"
}
}
}
}

