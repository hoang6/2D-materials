#! /usr/bin/perl

$nargs = @ARGV;
if($nargs != 5) {
    print "Usage: sampleRecords.pl <file> <^key> <start> <interval> <end>\n";
    exit(1);
}
$fname = $ARGV[0];
$key = $ARGV[1];
$a = $ARGV[2];
$d = $ARGV[3];
$z = $ARGV[4];

open(FIN, "<$fname") || die "ERROR: Cannot open $fname: $!\n";
$count = 0;
$flag = 0;
while(<FIN>) {
    chomp;
    if(/^$key/) {
	$count++;
	if($z >= 0) {
	    $flag = ($count >= $a && $count <= $z && ($count-$a)%$d == 0);
	}
	else {
	    $flag = ($count >= $a && ($count-$a)%$d == 0);
	}
    }
    print "$_\n" if $flag;
}
close(FIN);

