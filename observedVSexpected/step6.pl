#!/home/genkvg/perl/bin/perl -w
# made by Konstantin V. Gunbin aka GenKVG
# required OS: Linux
# required tools (installed in /usr/local/bin by default): no requirenments
# required script location: the directory containing rnd.* [see step4.pl], 'patterns' and step5_6.patts files
# run example: ./step6.pl EXPECTED step5_6.patts 

%IU=(
"A"=>"A",
"T"=>"T",
"G"=>"G",
"C"=>"C",
"R"=>"A|G",
"Y"=>"C|T",
"B"=>"C|G|T",
"D"=>"A|G|T",
"H"=>"A|C|T",
"V"=>"A|C|G",
"N"=>"A|T|G|C");

%comp=(
"-"=>"-",
"A"=>"T",
"T"=>"A",
"C"=>"G",
"G"=>"C",
"Y"=>"R",
"R"=>"Y",
"N"=>"N");

%PATTERNS=();
$lN=0;
open (P, "$ARGV[1]");
while (<P>)
{
chomp;
$lN++;
@P=split(/\s+/);
@patt=("$P[1]-$P[3]","$P[2]-$P[4]");
@patt0=split(/-/, $patt[0]);
@patt1=split(/-/, $patt[1]);
if ($P[0] eq "revcom")
{
foreach $p (@patt)
{
@p=split(//, $p);
$c=""; for ($q=$#p;$q>-1;$q--) {$c.=$comp{$p[$q]}}
push @cpatt, $c;
}
@cpatt0=split(/-/, $cpatt[0]);
@cpatt1=split(/-/, $cpatt[1]);
$PATTERNS{"$patt0[0]\t$patt0[1]\t$lN"}="$patt1[0]\t$patt1[1]";
$PATTERNS{"$cpatt0[0]\t$cpatt0[1]\t$lN"}="$cpatt1[0]\t$cpatt1[1]";
}
elsif ($P[0] eq "dir")
{
$PATTERNS{"$patt0[0]\t$patt0[1]\t$lN"}="$patt1[0]\t$patt1[1]";
}
}
close P;

%patternsKnow=();
open (O, "patterns");
while (<O>)
{
chomp;
@X=split(/\t/);
@x0=split(//, $X[0]);
@x1=split(//, $X[1]);
@x2=split(//, $X[2]);
@x3=split(//, $X[3]);

foreach $patt (sort {$a cmp $b} keys %PATTERNS)
{
@inpatt=split(/\t/, $patt);
@outpatt=split(/\t/, $PATTERNS{$patt});

@inpatt0=split(//, $inpatt[0]);
@inpatt1=split(//, $inpatt[1]);
@outpatt0=split(//, $outpatt[0]);
@outpatt1=split(//, $outpatt[1]);

$t0=0; for ($i=0;$i<3;$i++) {if ($IU{$inpatt0[$i]}=~/$x0[$i]/) {$t0+=1}}
$t1=0; for ($i=0;$i<3;$i++) {if ($IU{$outpatt0[$i]}=~/$x1[$i]/) {$t1+=1}}
$t2=0; for ($i=0;$i<3;$i++) {if ($IU{$inpatt1[$i]}=~/$x2[$i]/) {$t2+=1}}
$t3=0; for ($i=0;$i<3;$i++) {if ($IU{$outpatt1[$i]}=~/$x3[$i]/) {$t3+=1}}

if ($t0==3 and $t1==3 and $t2==3 and $t3==3)
    {
    $patternsKnow{"$X[0]\t$X[1]\t$X[2]\t$X[3]"}=$X[4];
    print STDERR "$X[0]\t$X[1]\t$X[2]\t$X[3] - $patt=>$PATTERNS{$patt}\n";
    }
}
}
close O;

###########
@dirs=(".");
###########
%FILES=();
$fN=0;
foreach $dir (@dirs)
{
open (L, "find $dir -type f -name 'rnd.*' |");
while (<L>)
{
chomp;
$file=$_;
@file=split(/\./, $file);
$fN++;
print STDERR "reading...$fN... $file $file[$#file] \n";
open (F, "$file");
while (<F>)
{
chomp;
@F=split(/\t/);
$patt="$F[0]\t$F[1]\t$F[2]\t$F[3]";
$numb=$F[4];
for ($i=0;$i<5;$i++) {shift @F}
if (exists $patternsKnow{$patt} and $numb==$patternsKnow{$patt})
{
foreach ($i=0;$i<$#F+1;$i++) {$FILES{$i}{$file[$#file]}+=$F[$i]}
}
}
close F;
}
close L;
}

open (M, "+>>$ARGV[0].mean.tsv");
open (S, "+>>$ARGV[0].sd.tsv");
$Mstr="$ARGV[1]";
$Sstr="$ARGV[1]";
foreach $pos (sort {$a<=>$b} keys %FILES)
{
@vals=();
foreach $it (keys %{$FILES{$pos}})
{
push @vals, $FILES{$pos}{$it};
}
print STDERR "$pos: @vals\n";
$mean=0; $N=0; foreach $v (@vals) {$mean+=$v; $N++} $mean=$mean/$N;
$sd=0; foreach $v (@vals) {$sd+=($v-$mean)*($v-$mean)} $sd=sqrt($sd/$N);
$Mstr.="\t$mean";
$Sstr.="\t$sd";
}
$Mstr.="\n";
$Sstr.="\n";
print M "$Mstr";
print S "$Sstr";
close M;
close S;
