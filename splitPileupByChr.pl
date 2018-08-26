#! /usr/bin/perl -w

use strict;
use Cwd;
use warnings;
use Getopt::Long;

sub usage {
        print <<USAGE;
usage:
        perl $0 [options]
        example: perl $0 
description:
        for splitting pileup file by chromosome
options:
		-help           : print help info;
		-in(str)        : input pileup file;
		-l(str)         : input chr list to split by;
		-out(str)       : output prefix;
		-gz(int)        : flag for output in gz type (*.gz), default: 0 (no);
USAGE
}

my ($flist, $glist, $gz_flag, $prefix, $chr_list, $out, $in, $mode, $help);
$gz_flag = 0;

GetOptions(
  "help"=>\$help,
  "in=s"=>\$in,
  "out=s"=>\$prefix,
  "l=s"=>\$chr_list,
  "gz=s"=>\$gz_flag,
);

&usage && exit if (defined $help || !defined $in || !defined $prefix || !defined $chr_list);

&showLog("start");
&showLog("read chr list $chr_list");
open(LIST, $chr_list) or die "can't open $chr_list";
my %hash_chr = ();
while (<LIST>) {
	chomp;
	next if /^#/;
	next if $_ eq "";
	my @line = split /\s+/;
	my $chr = $line[0];
	$hash_chr{$chr} = 0;
}

&showLog("read $in");
if ($in =~ /\.gz$/) {
	open(IN, "gzip -dc $in |") or die "can't open file $in";
} else {
	open(IN, "$in") or die "can't open file $in";
}

my $STAND_LENGTH = 0;

my %hash_filehandle = ();
while (<IN>) {
	chomp;
	next if ($_ eq "");
	next if /^#/;
	my @line = split /\s+/;
	my $chr = $line[0];
	next if (scalar(@line) < 6);
	next if (!exists($hash_chr{$chr}));
	if ($STAND_LENGTH == 0) {
		$STAND_LENGTH = scalar(@line);
	}
	next if (scalar(@line) != $STAND_LENGTH && $STAND_LENGTH != 0);
	if (!exists($hash_filehandle{$chr})) {
		&showLog("output to $chr");
		if ($gz_flag) {
			open($hash_filehandle{$chr}, "|gzip -c >$prefix.$chr.pileup.gz") or die "can't write file $prefix.$chr.pileup.gz";
		} else {
			open($hash_filehandle{$chr}, ">$prefix.$chr.pileup") or die "can't write file $prefix.$chr.pileup";
		}
	}
	print {$hash_filehandle{$chr}} "$_\n";
}
close IN;
&showLog("all done");

sub showLog {
	my ($info) = @_;
	my @times = localtime; # sec, min, hour, day, month, year
	print STDERR sprintf("[%d-%02d-%02d %02d:%02d:%02d] %s\n", $times[5] + 1900, $times[4] + 1, $times[3], $times[2], $times[1], $times[0], $info);
}

