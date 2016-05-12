#! /usr/bin/perl

## DEPENDENCY:
 # "./getGPHDynMat"
 # "./makeGPHGroundState"
## REMARK:
 # This script executes getGPHDynMat for graphenes with 
 # different lattice parameter a (in equilibrium a = 2.46 A)

###############################################################
############################ ARGS #############################
###############################################################
$nargs = $#ARGV+1;
if($nargs != 11) {
    print "Usage: batch.getGPHDynMat.pl <a0> <da> <na> <nchunk> <ncell1D> <nq1D[0]> <nq1D[1]> <work_dir> <file_mol> <file_dynMat> <test>\n";
    exit(1);
}
$a0 = $ARGV[0];
$da = $ARGV[1];
$na = $ARGV[2];
$nchunk = $ARGV[3];
$ncell1D = $ARGV[4];
@nq1D = ($ARGV[5],$ARGV[6]);
$workDir = $ARGV[7];
$molFile = $ARGV[8];
$dynMatFile = $ARGV[9];
$test = $ARGV[10];
$tag = "temp";

###############################################################
############################ MAIN #############################
###############################################################
# main loop
$cmd = "";
$format = "%0".length($nchunk)."s";
$chunkSize = int($na/$nchunk);
$chunkSize++ if ($chunkSize < $na/$nchunk); 
foreach $i (0..($nchunk-1)) {
    $tail = sprintf($format,$i);
    $myWorkDir = ($workDir).".".($tail);
    $myMolFile = ($myWorkDir)."/".($molFile);
    $myDynMatFile = ($myWorkDir)."/".($dynMatFile);
    $myTmpMolFile = ($myMolFile).".".($tag);
    $myTmpDynMatFile = ($myDynMatFile).".".($tag);
    $myJunkFile = ($myWorkDir)."/junk.".($tag);
    $cmd .= "mkdir -p $myWorkDir; ";
    $jmin = $chunkSize*$i;
    $jmax = $chunkSize*($i+1)-1;
    $jmax = $na-1 if ($i == ($nchunk-1));
    foreach $j ($jmin..$jmax) {
	$a = $a0+$j*$da;
	$cmd .= "./makeGPHGroundState $a $ncell1D 1 $myTmpMolFile > $myJunkFile; ";
	$cmd .= "./getGPHDynMat $myTmpMolFile 1.0e-5 1.0e-6 1 $nq1D[0] $nq1D[1] 0 > $myTmpDynMatFile; ";
	$cmd .= "cat $myTmpMolFile >> $myMolFile; ";
	$cmd .= "cat $myTmpDynMatFile >> $myDynMatFile; ";
    }
    $cmd .= "rm -f $myTmpMolFile $myTmpDynMatFile $myJunkFile\n\n";
}
# end main loop

print "$cmd";
system($cmd) if (!$test);
