#!/usr/bin/perl -X
# made by Konstantin V. Gunbin aka GenKVG
# required OS: Linux
# required tools (installed in /usr/local/bin by default): no requirenments
# required script location: the directory containing OUTDISTANCES.random.* [see step3.pl] and OUTDISTANCES.real.* [see step4.pl] files
# run example: ./step5.pl OUTDISTANCES.real.* OUTDISTANCES.random.* distance_to_burnin[e.g. 2] > RESULTS.stats 

%REAL=();
$l=0;
open (I, "$ARGV[0]");
while (<I>)
{
$l++;
chomp;
@I=split(/\t/);
if ($l==1)
{
shift @I; shift @I; shift @I;
@heads=@I;
}
else
{
#SampleID        Dinucleotide    Change
shift @I;
$dn=shift @I;
$chg=shift @I;
$str="$chg:$dn";
for ($i=0;$i<$#I+1;$i++)
    {
    if ($heads[$i]>=$ARGV[2] and $I[$i]>0) {for ($j=1;$j<$I[$i]+1;$j++) {push @{$REAL{$str}{3}}, $heads[$i]; # print STDERR "@{$REAL{$str}{3}}\n"
    }}
    if ($heads[$i]<=-$ARGV[2] and $I[$i]>0) {for ($j=1;$j<$I[$i]+1;$j++) {push @{$REAL{$str}{5}}, -($heads[$i]); #print STDERR "@{$REAL{$str}{5}}\n"
    }}
    }
}
}
close I;

%RANDOM=();
$l=0;
open (I, "$ARGV[1]");
while (<I>)
{
$l++;
chomp;
@I=split(/\t/);
if ($l==1)
{
shift @I; shift @I;
@heads=@I;
}
else
{
#Dinucleotide    Change
$dn=shift @I;
$chg=shift @I;
$str="$chg:$dn";
for ($i=0;$i<$#I+1;$i++)
    {
    if ($heads[$i]>=$ARGV[2] and $I[$i]>0) {for ($j=1;$j<$I[$i]+1;$j++) {push @{$RANDOM{$str}{3}}, $heads[$i]; #print STDERR "@{$RANDOM{$str}{3}}\n"
    }}
    if ($heads[$i]<=-$ARGV[2] and $I[$i]>0) {for ($j=1;$j<$I[$i]+1;$j++) {push @{$RANDOM{$str}{5}}, -($heads[$i]); #print STDERR "@{$RANDOM{$str}{5}}\n"
    }}
    }
}
}
close I;

foreach $str (sort {$a cmp $b} keys %RANDOM)
{
@str=split(/:/, $str);
$L="$str[0]\t$str[1]";
@vals=@{$REAL{$str}{5}}; #print STDERR "$L\treal5\t@vals\n";
&statscount(); $m5real=$mean; $gm5real=$gm; $hm5real=$hm; $md5real=$median; $Q5_25real=$q25; $Q5_75real=$q75; $n5real=$n+1;
@vals=@{$REAL{$str}{3}}; #print STDERR "$L\treal3\t@vals\n";
&statscount(); $m3real=$mean; $gm3real=$gm; $hm3real=$hm; $md3real=$median; $Q3_25real=$q25; $Q3_75real=$q75; $n3real=$n+1;
@vals=@{$RANDOM{$str}{5}}; #print STDERR "$L\trand5\t@vals\n";
&statscount(); $m5rand=$mean; $gm5rand=$gm; $hm5rand=$hm; $md5rand=$median; $Q5_25rand=$q25; $Q5_75rand=$q75; $n5rand=$n+1;
@vals=@{$RANDOM{$str}{3}}; #print STDERR "$L\trand3\t@vals\n";
&statscount(); $m3rand=$mean; $gm3rand=$gm; $hm3rand=$hm; $md3rand=$median; $Q3_25rand=$q25; $Q3_75rand=$q75; $n3rand=$n+1;
$L.="\t$n5real\t$m5real\t$gm5real\t$hm5real\t$md5real\t$Q5_25real\t$Q5_75real\t$n5rand\t$m5rand\t$gm5rand\t$hm5rand\t$md5rand\t$Q5_25rand\t$Q5_75rand\t$n3real\t$m3real\t$gm3real\t$hm3real\t$md3real\t$Q3_25real\t$Q3_75real\t$n3rand\t$m3rand\t$gm3rand\t$hm3rand\t$md3rand\t$Q3_25rand\t$Q3_75rand";
print "$L\n";
}

sub statscount
{
if ($#vals==-1) {$mean=$hm=$gm=$median=$q25=$q75="NA"; $n=0}
else
{
$mean=0;
$hm=0;
$n=$#vals;
foreach $v (@vals) {$mean+=$v} $mean=$mean/($n+1);
foreach $v (@vals) {$hm+=(1/$v)} $hm=($n+1)/$hm;
$gm=sqrt($hm*$mean);
@svals=sort {$a<=>$b} @vals;
$median=$svals[int($n/2)];
$q25=$svals[int($n/4)];
$q75=$svals[int((3*$n)/4)];
}
}
