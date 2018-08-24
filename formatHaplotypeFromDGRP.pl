#! /usr/bin/perl -w

use strict;
use Cwd;
use warnings;
use Getopt::Long;

sub usage {
    print <<USAGE;
usage:
    perl $0 [options]
    example: perl $0 -help
options:
    -help       : print help info;
    -i(str)     : input file containing genotype from DGRP2;
    -n(int)     : number of samples to generate, default 100;
    -m(int)     : number of SNPs to generate, randomly pick up from the input file (-i), default 100000;
    -o(str)     : output file;
USAGE
}

my ($help, $in, $n, $m, $out);
$n = 100;
$m = 100000;

GetOptions(
"help"=>\$help,
"i=s"=>\$in,
"o=s"=>\$out,
"n=s"=>\$n,
"m=s"=>\$m,
);

if (defined $help || !defined $in || !defined $out) {
    &usage();
    exit;
}

sub showLog {
    my ($info) = @_;
    my @times = localtime; # sec, min, hour, day, month, year
    print STDERR sprintf("[%d-%02d-%02d %02d:%02d:%02d] %s\n", $times[5] + 1900, $times[4] + 1, $times[3], $times[2], $times[1], $times[0], $info);
}

&showLog("start");

&showLog("read input file $in");
my $snpCount = `grep SNP $in | wc -l`;
my $colCount = `head -n 1 $in | awk '{print NF}'`;
my $offset = 9;
my $sampleNumber = $colCount - $offset;
chomp($snpCount);
open(IN, "$in") or die "can't open file $in";
open(OUT, ">$out") or die "can't write to file $out";
my $count = 0;
srand(time);
&showLog("Total SNPs: $snpCount; Total samples: $sampleNumber");
while (<IN>) {
    chomp;
    next if ($_ eq "");
    next if /^#/;
    $count++;
    &showLog("processed $count lines") if ($count % 10000 == 0);
    my $line = $_;
    my @lines = split /\s/, $line;
    my ($chr, $pos, $id, $ref, $alt, $refc, $altc, $qual) = @lines;
    next if (!($id =~ /SNP/));
    my $output;
    my ($refcount, $altcount) = (0, 0);
    if (rand(1) <= $m/($snpCount + 1)) {
        for (my $i = 0; $i < $n; $i++) {
            my $index = ($i % $sampleNumber) + $offset; # if the number of samples exceeds the available lines from DGRP, reuse the lines
            if ($lines[$index] eq "0") {
                $output .= "$ref$ref ";
                $refcount++;
            } elsif ($lines[$index] eq "2") {
                $output .= "$alt$alt ";
                $altcount++;
            } elsif ($lines[$index] eq "-") {
                if (rand(1) <= 0.5) {
                    $output .= "$ref$ref ";
                    $refcount++;
                } else {
                    $output .= "$alt$alt ";
                    $altcount++;
                }
            }
        }
        $output =~ s/\s+$//;
        if ($refcount < $altcount) {
            print OUT join("\t", $chr, $pos, $ref, "$alt/$ref", $output), "\n";
        } else {
            print OUT join("\t", $chr, $pos, $ref, "$ref/$alt", $output), "\n";
        }
    }
}
close IN;
close OUT;
&showLog("read and output finished");

### output haplotype format:
#2L      686891    G      G/A    AG GG GG GG GG
#2L      936681    A      A/G    GG AA AA AA AA
#2L      861026    T      A/T    TT AT AA AA TT
#2L      966618    C      T/C    TC TC TT TT TT
#2R      134298    A      A/C    AC AC CC CC CC
#The file consists of exactly five columns which are separated by a 'tab'!

#col1: chromosome
#col2: position
#col3: reference character
#col4: allelels of the SNP, in the foramt "major-allele/minor-allele"
#col5: a space separated list of genotypes. The first genotype refers to the first specimen, the second to the second specimen and so on..
