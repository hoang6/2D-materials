#! /usr/bin/perl

$wdir = "CNTs";
open(FMAN,">job.man");
open(FOPT,">job.opt");
print FMAN "mkdir $wdir\n";
foreach $m (3..12) {
    $n = int(sqrt(3.0)*$m+0.5);
    $tag = sprintf("%02d",($m-2));
    $arg = "$m $m -$n $n"; $fname = "$wdir/aCNT$tag";
    print FMAN "./manufacture -gm CNT \"$arg\" > $fname\n";
    print FOPT "./testOptim REBO $fname $fname.relaxed > junk.aCNT$tag\n\n";
    $arg = "$n 0 -$m ".2*$m; $fname = "$wdir/zCNT$tag";
    print FMAN "./manufacture -gm CNT \"$arg\" > $fname\n";
    print FOPT "./testOptim REBO $fname $fname.relaxed > junk.zCNT$tag\n\n";
}
close(FOPT);
close(FMAN);
