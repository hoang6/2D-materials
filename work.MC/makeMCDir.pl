#! /usr/bin/perl

## DEPENDENCY:
 # "../clt"
 # "../work.util/simpleMPIDriver"
 # "../work.util/manufacture"
 # "../work.util/clone.pl"
 # "./runMC<ensemble>"
 # "./analMC<ensemble>"
 # "./analGPH"
 # "./analCNT"
 # "./repairMCLogs.pl"
 # "./resumeMC.pl"
 # "./bat.MC<ensemble>.pl"
 # "./in.bat.MC<ensemble>.sample"
 # "./nRecords.pl"
 # "./sampleRecords.pl"
## REMARK:
 # This script makes a MC simulation directory,
 # i.e. copy everything needed to this directory.

###############################################################
############################ MAIN #############################
###############################################################
# read args
$nargs = $#ARGV+1;
if($nargs != 2 && $nargs != 3) {
    print "Usage: makeMCDir <ensemble_name> <dir_name> [test]\n";
    exit 1;
}
$ensemble = $ARGV[0];
$dir = $ARGV[1];
if($nargs == 3) {
    $test = $ARGV[2];
}
else {# default
    $test = 0;
}

# make dir
$cmd = "";
$cmd .= "mkdir -p $dir\n";
$cmd .= "cp ../clt $dir\n";
$cmd .= "cp ../work.util/simpleMPIDriver $dir\n";
$cmd .= "cp ../work.util/manufacture $dir\n";
$cmd .= "cp ../work.util/clone.pl $dir\n";
$cmd .= "cp ./runMC".($ensemble)." $dir\n";
$cmd .= "cp ./analMC".($ensemble)." $dir\n";
$cmd .= "cp ./analGPH $dir\n";
$cmd .= "cp ./analCNT $dir\n";
$cmd .= "cp ./repairMCLogs.pl $dir\n";
$cmd .= "cp ./resumeMC.pl $dir\n";
$cmd .= "cp ./bat.MC".($ensemble).".pl $dir\n";
$cmd .= "cp ./in.bat.MC".($ensemble).".sample $dir\n";
$cmd .= "cp ./nRecords.pl $dir\n";
$cmd .= "cp ./sampleRecords.pl $dir\n";

print "Execute\n$cmd";
system($cmd) if !$test;
