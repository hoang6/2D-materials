#! /usr/bin/perl

## customize here
$m = 12;    # GPH_size = (m,m,-h,h). Unit: 2.46A
$nL = 1;   # NO. of layers
$dz = 0.5; # indenter step size. Unit: A
$nS = 20;  # NO. of steps 

## parameters
$h = int(sqrt(3)*$m);            # GPH_size = (m,m,-h,h)
$r = int($m*2.46*sqrt(3)/2.0)-2; # hole radius

## generate cmd
$mfile = "GPH${nL}L";
$hfile = "${mfile}.Hole";
$cmd = "";
$cmd .= "./manufacture -am GPH '$m $m -$h $h' > out.log\n";
$cmd .= "./testGPHnLa GPH* $mfile $nL >> out.log\n";
$cmd .= "./testGPHonHole $mfile $hfile $r 0.0 >> out.log\n\n";
#print("$cmd");
system("$cmd");

open(FIN,"<out.log");
while(<FIN>) {
    chomp;
    if(/atomID = ([0-9]+)/) {
	$atomID = $1; last;
    }
}
close(FIN);
#print("atomID = $atomID\n");

$cmd = "";
for($i = 0; $i <= $nS; $i++) {
    $cfile = "${hfile}.d".sprintf("%03d",$i);
    if($i == 0) {
	$cmd .= "cp $hfile $cfile; ";
    }
    else {
	$cmd .= "./testDispAtom $pfile $cfile $atomID 0.0 0.0 -$dz; "
    }
    $pfile = $cfile;
    $cmd .= "./testOptim LCBOPII $cfile $cfile >> out.log; ";
    $cmd .= "rm *.pro *.mes; \n\n";
}
print($cmd);
