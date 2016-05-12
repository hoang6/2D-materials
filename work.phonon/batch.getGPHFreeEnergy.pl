#! /usr/bin/perl

## DEPENDENCY
 # "./getGPHFreeEnergy"
## REMARK
 # This script executes getGPHFreeEnergy for outputs of batch.getGPHDynMat.pl

###############################################################
############################ ARGS #############################
###############################################################
$nargs = $#ARGV+1;
if($nargs != 9) {
    print "Usage: batch.getGPHFreeEnergy.pl <T0> <dT> <nT> <nchunk> <work_dir> <file_mol> <file_dynMat> <file_aT> <test>\n";
    exit(1);
}
$T0 = $ARGV[0];
$dT = $ARGV[1];
$nT = $ARGV[2];
$nchunk = $ARGV[3];
$workDir = $ARGV[4];
$molFile = $ARGV[5];
$dynMatFile = $ARGV[6];
$aTFile = $ARGV[7];
$test = $ARGV[8];

###############################################################
############################ MAIN #############################
###############################################################
# main loop
$cmd = "";
$format = "%0".length($nchunk)."s";
foreach $i (0..($nchunk-1)) {
    $tail = sprintf($format,$i);
    $myWorkDir = ($workDir).".".($tail);
    $myMolFile = ($myWorkDir)."/".($molFile);
    $myDynMatFile = ($myWorkDir)."/".($dynMatFile);
    $myATFile = ($myWorkDir)."/".($aTFile);
    $cmd .= "./getGPHFreeEnergy $myMolFile $myDynMatFile $T0 $dT $nT > $myATFile\n\n";
}

print "$cmd";
system($cmd) if (!$test);
