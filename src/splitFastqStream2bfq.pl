#!/usr/bin/perl
use strict;
use warnings;

my $usage = "USAGE: splitFastqStream2bfq.pl <size> <fastq_prefix> \[<extension> defaults to bfq\]\n";
my $count = shift or die $usage;
my $prefix = shift or die $usage;
my $ext = shift;
$ext = 'bfq' unless($ext);

my $fc = 1;
my $c = 0;
my $one_line_seen;

open(my $out_h, '>', "${prefix}_${fc}.${ext}") or die $!;
while (my $line = <>) {
    $one_line_seen++;
    if ($c == $count) {
        close $out_h;
        $fc++;
        open($out_h, '>', "${prefix}_${fc}.${ext}") or die $!;
        $c = 0;
    }
    print $out_h $line;
    $line = <>;
    print $out_h $line;
    $line = <>;
    print $out_h $line;
    $line = <>;
    print $out_h $line;
    $c++;
}
close($out_h);
if ($one_line_seen) {
    print $fc;
}
else {
    unlink "${prefix}_${fc}.${ext}";
    print 0;
}
exit;
