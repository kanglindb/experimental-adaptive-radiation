#! /usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use POSIX;

sub usage {
    print <<USAGE;
usage:
    perl $0 [options]
    example: perl $0 -help
description:
    function for picking p-value threshold of CMH test for each category based on starting minor allele frequencies
    (2018-08-15) first version
options:
    -help       : print help info;
    -i(str)     : input file CMH test result (output from popoolation2/cmh-test.pl);
    -fdr(int)   : coresponding false dicovery rate, default 1;
USAGE
}

my ($help, $in, $fdr);
$fdr = 1;

GetOptions(
"help"=>\$help,
"i=s"=>\$in,
"fdr=s"=>\$fdr,
);

if (defined $help || !defined $in) {
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

open(IN, "$in") or die "can't open file $in";
my $count = 0;
my (@countGroup) = ();
my (%hash_pvalue) = ();
&showLog("read file $in");
while (<IN>) {
    chomp;
    next if ($_ eq "");
    next if /^#/;
    $count++;
    &showLog("processed $count lines") if ($count % 100000 == 0);
    my $line = $_;
    my @lines = split /\s/, $line;
    my ($chr, $pos, $ref, @sync) = @lines;
    my $freq = &getFreq($sync[0]);
    my $index = floor($freq/10);
    $index = 4 if ($index > 4);
    $hash_pvalue{$index}{$sync[-1]} = 0;
    $countGroup[$index]++;
}
close IN;

foreach my $group (sort {$a <=> $b} keys %hash_pvalue) {
    print "Category ", $group * 10, "-", ($group + 1) * 10, ": ";
    my @tmpArray = sort {$a <=> $b} keys %{$hash_pvalue{$group}};
    my $index = int($fdr/100 * $countGroup[$group]);
    $index = 1 if ($index < 1);
    print $tmpArray[$index - 1], "\n";
}

&showLog("all done");

sub getFreq {
    my $str = shift;
    # sync format: allele frequencies for all populations in the form A-count:T-count:C-count:G-count:N-count:deletion-count
    my ($na, $nt, $nc, $ng) = split(":", $str);
    my @tmpArray = ($na, $nt, $nc, $ng);
    my $total = $na + $nt + $nc + $ng;
    if ($total == 0) {
        print STDERR "WARNING, parsing SYNC file format failed!\n";
        return 0;
    }
    @tmpArray = sort {$b <=> $a} @tmpArray;
    return int($tmpArray[1]/$total * 100);
}

