#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

sub usage {
        print <<USAGE;
usage:
        perl $0 [options]
        example: perl $0 in1.file in2.file in3.file ... > merge.file
description:
        function for merging files by keys.
        (2017-12-20) first version
options:
        -help      : print help info;
        -k(str)    : key index default, multiple indices separated by COMMA (e.g 1,2), default: 1;
USAGE
}

my $keyPos = 1;
my $help;
my $repeatFlag = 0; ### flag for outputing lines shared by all input files only
my $flag = 0; ### flag for outputing the original content from each input file, by default, the index column(s) will be printed only once
my $outPos; ### if defined, only output the #$outPos columns from each file 

GetOptions(
    "help"=>\$help,
    "k=s"=>\$keyPos,
    "l=s"=>\$outPos,
    "d=s"=>\$repeatFlag,
    "f=s"=>\$flag,
);

if ($#ARGV < 0 || defined $help || !defined $keyPos) {
    usage();
    exit;
}

my %hash;
my %hash_count;
my $i = 0;
my $count = 0;
my @files = @ARGV;
foreach my $file (@files) {
    &merge($file, $count++) if (-e $file); 
}

&showLog("output file");
print STDERR "count: $count\n";
foreach my $key (sort {$a cmp $b} keys %hash) {
    next if ($repeatFlag && $hash_count{$key} <= $#ARGV);
    if ($flag) {
        for (my $j = 0; $j <= $#ARGV - 1; $j++) {
            print $hash{$key}[$j], "\t";
        }
        print $hash{$key}[$#ARGV], "\n";
    } else {
        print $key;
        for (my $j = 0; $j <= $#ARGV; $j++) {
            print "\t", $hash{$key}[$j];
        }
        print "\n";
    }
}
&showLog("done");

sub merge {
    my $file = shift;
    my $index = shift;
    &showLog("read file $file");
    open FILE, $file or print "cant open file $file!\n";
    while (<FILE>) {
    	chomp;
		my @line = split("\t", $_);
		my $key = &getKey($keyPos, \@line);
		my $supContent = "-";
		for (my $i = 2; $i < scalar(@line); $i++) {
			$supContent .= "\t-";
		}
		my $content = join("\t", @line);
		if (!($flag)) {
			my @tmp_line = @line;
			my @keys = split(",", $keyPos);
			foreach my $index (@keys) {
				$tmp_line[$index - 1] = "";
			}
			$content = join("\t", @tmp_line);
		}
		if (defined $outPos) {
			$content = &getKey($outPos, \@line);
			my @tmp_line = split("", "----------------------------------------------------------------------");
			$supContent = &getKey($outPos, \@tmp_line);
		}
		if (exists($hash{$key})) {
			$hash_count{$key}++;
			$hash{$key}[$index] = $content;
		} else {
			for (my $j = 0; $j <= $#ARGV; $j++) {
				push(@{$hash{$key}}, $supContent);
			}
			$hash_count{$key} = 1;
			$hash{$key}[$index] = $content;
		}
	}
	close FILE;
}

sub getKey {
    my ($keyStr, $arry_ref) = @_;
    my @keys = split(",", $keyStr);
    my $tmp_key = "";
    my @tmp;
    foreach my $index (@keys) {
        push(@tmp, $arry_ref->[$index - 1]);
    }
    $tmp_key = join("\t", @tmp);
    return $tmp_key;
}

sub showLog {
	my ($info) = @_;
	my @times = localtime; # sec, min, hour, day, month, year
	print STDERR sprintf("[%d-%02d-%02d %02d:%02d:%02d] %s\n", $times[5] + 1900, $times[4] + 1, $times[3], $times[2], $times[1], $times[0], $info);
}

