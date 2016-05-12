#! /usr/bin/perl

$nargs = @ARGV;
if($nargs != 2) {
    print "Usage: nRecords.pl <file> <^key>\n";
    exit(1);
}
$fname = $ARGV[0];
$key = $ARGV[1];

open(FIN, "<$fname") || die "ERROR: Cannot open $fname: $!\n";
$count = 0;
while(<FIN>) {
    chomp;
    $count++ if(/^$key/);
}
close(FIN);

print "$count\n";
