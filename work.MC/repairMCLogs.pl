#! /usr/bin/perl

$nargs = $#ARGV+1;
if($nargs < 7 || $nargs > 9 ) {
    print "Usage: repairMCLogs.pl <in.log> <action.log> <state.log> <stat.log> <info.log> <max # of mols)> <sample_freq> [digits] [tag]\n";
    exit(1);
}
$inFile = $ARGV[0];
$actionFile = $ARGV[1];
$stateFile = $ARGV[2];
$statFile = $ARGV[3];
$infoFile = $ARGV[4];
$nmols = $ARGV[5];
$freq = $ARGV[6];
if($nargs > 7) {
    $digits = $ARGV[7];
}
else {
    $digits = -1;
}
if($nargs > 8) {
    $tag = $ARGV[8];
}
else {
    $tag = ".temp";
}

$inFileTemp = ($inFile).($tag);
$actionFileTemp = ($actionFile).($tag);
$stateFileTemp = ($stateFile).($tag);
$statFileTemp = ($statFile).($tag);
$infoFileTemp = ($infoFile).($tag);

$cmdMV = "mv $inFileTemp $inFile; mv $actionFileTemp $actionFile; mv $stateFileTemp $stateFile; mv $statFileTemp $statFile; mv $infoFileTemp $infoFile;";
$cmdRM = "rm $inFileTemp $actionFileTemp $stateFileTemp $statFileTemp $infoFileTemp;";

system("touch $inFileTemp $actionFileTemp $stateFileTemp $statFileTemp $infoFileTemp;");

## main
# repair stateFile
$indRange = ($freq-1).":".($freq).":".($nmols-1);
$digitsArg = "";
if($digits >= 0) { $digitsArg = "-d $digits"; }
$result = `./extract -k $stateFile [$indRange] $digitsArg -o $stateFileTemp`;
$result =~ m/Output\s+([0-9]+)\s+molecules/i;
$nmols = $1;
if($nmols == 0) {
    print "ERROR: Cannot have zero molecules\n";
    system($cmdRM);
    exit(2);
}

# repair actionFile, statFile, infoFile
@arrayFile = ($actionFile, $statFile, $infoFile);
@arrayFileTemp = ($actionFileTemp, $statFileTemp, $infoFileTemp);
foreach $k (0,1,2) {
    $aFile = $arrayFile[$k];
    $tFile = $arrayFileTemp[$k];
    open(FOUT, ">$tFile") || die "ERROR: Cannot open $tFile: $!\n";
    open(FIN, "<$aFile")  || die "ERROR: Cannot open $aFile: $!\n";
    $count = 0;
    $flag = 0;
    while($tline = <FIN>) {
	chomp $tline;
	if($tline eq "") {
	    if(($count+1)%($freq) == 0) { 
		print FOUT "$tline\n";
		if(($count+1)/($freq) == $nmols) { $flag = 1; last; }
	    }
	    $count++;
	}
	else {
	    print FOUT "$tline\n";
	}
    }
    close(FIN)  || die "ERROR: Cannot close $aFile: $!\n";
    close(FOUT) || die "ERROR: Cannot close $tFile: $!\n";
    if(!$flag) {
	print "ERROR: $aFile has < $nmols blank lines\n";
	system($cmdRM);
	exit(3);
    }
}

# repair inFile
open(FIN, "<$inFile") || die "ERROR: Cannot open $inFile: $!\n";
while($tline = <FIN>) {
    if($tline =~ m/^dStateStep\s*([0-9]+)/i) {
	$dStateStep = $1;
    }
}
$dStateStep = $dStateStep*$freq;
$eStep = $nmols*$dStateStep;
close(FIN)  || die "Cannot close $inFile: $!\n";
open(FOUT, ">$inFileTemp") || die "ERROR: Cannot open $inFileTemp: $!\n";
open(FIN, "<$inFile")      || die "ERROR: Cannot open $inFile: $!\n";
while($tline = <FIN>) {
    if($tline =~ m/^(eStep)(\s+)/i) {
	print FOUT "$1$2$eStep\n";
    }
    elsif($tline =~ m/^(dStateStep)(\s+)/i) {
	print FOUT "$1$2$dStateStep\n";
    }
    else {
	print FOUT $tline;
    }
}
close(FIN)  || die "Cannot close $inFile: $!\n";
close(FOUT) || die "Cannot close $inFileTemp: $!\n";

system($cmdMV);
