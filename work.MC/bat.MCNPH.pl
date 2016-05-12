#! /usr/bin/perl

## DEPENDENCY:
 # "./simpleMPIDriver"
 # "./runMCNPH"
 # "./extract"
 # "./repairMCLogs.pl"
 # "./resumeMC.pl"
 # "./analMCNPH"
 # "./analGPH or ./analCNT"
## REMARK:
 # This script sets up the NPH MC simulations

###############################################################
############################ ARGS #############################
###############################################################
$nargs = $#ARGV+1;
if($nargs != 2 && $nargs != 3) {
    print "Usage: bat.runMCNPH.pl <in.bat> <task(new/repair/resume/anal)> [tag]\n";
    print "       new      <init.mol>                    (> job.script)\n";
    print "       repair   <max (nmols)> <sample freq>   (> repair.script)\n";
    print "       resume   <rStep> <eStep>               (> resume.script)\n";
    print "       anal...  <skip(nmols)> <drec(nmols)>   (> anal.script)\n";
    exit(1);
}
$inBatFile = $ARGV[0];
$task = $ARGV[1];
if($nargs == 3) {
    $tag = $ARGV[2];
}
else {
    $tag = ".temp";
}

###############################################################
############################ SETUP ############################
###############################################################
# default setup
$simTag        = "MCNPH";
$potTag        = "REBO";
$agentTag      = "Default";
$rootDir       = ".";
@simDirArray   = ("exp1");
@caseDirArray  = ("case1");
@enthalpyArray = ("-6.6");
@pressureArray = ("0.0");
@extParaName   = ();
@extParaArray  = ();
$rStep         = 100;
$eStep         = 200;
$hStep         = 1;
$dAction       = 1;
$dStateStep    = 10;
$dStatStep     = 10;
$dInfoStep     = 10;
$inFile        = "in.log";
$actionFile    = "action.log";
$stateFile     = "state.log";
$statFile      = "stat.log";
$infoFile      = "info.log";
$junkFile      = "junk.log";
$digits        = 9;

# read the setup file
open(FINBAT, "<$inBatFile") || die "ERROR: Cannot open $inBatFile: $!\n";
while($tline = <FINBAT>) {
    chomp($tline);
    if($tline =~ m/^#/) {
	next;
    }
    elsif($tline =~ m/^simTag\s+(\S+)/i) {
	$simTag = $1;
    }
    elsif($tline =~ m/^potential\s+(\S+)/i) {
	$potTag = $1;
    }
    elsif($tline =~ m/^rStep\s+([0-9]+)/i) {
	$rStep = $1;
    }
    elsif($tline =~ m/^eStep\s+([0-9]+)/i) {
	$eStep = $1;
    }
    elsif($tline =~ m/^hStep\s+([0-9]+)/i) {
	$hStep = $1;
    }
    elsif($tline =~ m/^dAction\s+([0-9])/i) {
	$dAction = $1;
    }
    elsif($tline =~ m/^dStateStep\s+([0-9]+)/i) {
	$dStateStep = $1;
    }
    elsif($tline =~ m/^dStatStep\s+([0-9]+)/i) {
	$dStatStep = $1;
    }
    elsif($tline =~ m/^dInfoStep\s+([0-9]+)/i) {
	$dInfoStep = $1;
    }
    elsif($tline =~ m/^inFile\s+(\S+)/i) {
	$inFile = $1;
    }
    elsif($tline =~ m/^actionFile\s+(\S+)/i) {
	$actionFile = $1;
    }
    elsif($tline =~ m/^stateFile\s+(\S+)/i) {
	$stateFile = $1;
    }
    elsif($tline =~ m/^statFile\s+(\S+)/i) {
	$statFile = $1;
    }
    elsif($tline =~ m/^infoFile\s+(\S+)/i) {
	$infoFile = $1;
    }
    elsif($tline =~ m/^junkFile\s+(\S+)/i) {
	$junkFile = $1;
    }
    elsif($tline =~ m/^digits\s+([0-9]+)/i) {
	$digits = $1;
    }
    elsif($tline =~ m/^rootDir\s+(\S+)/i) {
	$rootDir = $1;
    }
    elsif($tline =~ m/^simDir\s+(\S.*)/i) {
	@simDirArray = split(/\s+/,$1);
    }
    elsif($tline =~ m/^tableBegin:\s*(\S+)/i) {
	$agentTag = $1;
	@caseDirArray = ();
	@enthalpyArray = ();
	@pressureArray = ();
	@extraParaArray = ();
	@extParaName = ();
	@extParaArray = ();
	$tline = <FINBAT>; # read the table head
	chomp($tline);
	@tableHead = split(/\s+/,$tline);
	while($tline = <FINBAT>) {
	    chomp($tline);
	    if($tline =~ m/^tableEnd/i) {
		last;
	    }
	    @fields = split(/\s+/,$tline);
	    if($#fields >= 2) {
		push(@caseDirArray,$fields[0]);
		push(@enthalpyArray,$fields[1]);
		push(@pressureArray,$fields[2]);
		$kmax = ($#fields-2);
		for($k = 1; $k <= $kmax; $k++) {
		    push(@extParaArray,$fields[$k+2]);
		}
	    }
	}#end while
	for($k = $kmax; $k >= 1; $k--) {
	    push(@extParaName,$tableHead[$#tableHead-($k-1)]);
	}
    }
}
close(FINBAT) || die "ERROR: Cannot close $inBatFile: $!\n";

###############################################################
############################ MAIN #############################
###############################################################
# glue dir names
@dirArray = ();
foreach $k (@simDirArray) {
    foreach $j (@caseDirArray) {
	push(@dirArray, ($rootDir)."/".($k)."/".($j));
    }
}

# new mode:
if($task =~ m/^new\s+(\S+)/i) {
    $mFile = $1;
    $nrep = $#caseDirArray+1;
    $count = 0;
    foreach $count (0..$#dirArray) {
	$cdir = $dirArray[$count];
	$seed = int(rand(1.0e9));
	$iFile = ($cdir)."/".($inFile);
	$jFile = ($cdir)."/".($junkFile);
	$pressure = $pressureArray[$count%$nrep];
	$enthalpy = $enthalpyArray[$count%$nrep];
	$approx_kT = ($enthalpy+7.3557)/1.5;
	$dxmax = 5.5*sqrt(-log(0.50)*$approx_kT)/31.0;
	$dhmax = 4.8*sqrt(-log(0.50)*$approx_kT)/20.0;
	system("mkdir -p ".($cdir));
	if(!(-d $cdir)) {
	    print "WARNING: Cannot make dir $cdir\n";
	    next;
	}
	open(FIN, ">$iFile") || die "ERROR: Cannot open $iFile: $!\n";
	print FIN "seed                $seed\n";
	print FIN "potential           $potTag\n";
	print FIN "MCAgent             $agentTag\n";
	print FIN "pressure            $pressure\n";
	print FIN "enthalpy            $enthalpy\n";
	$nName = $#extParaName+1;
	for($k = 0; $k < $nName; $k++) {
	    $tname = sprintf("%-20s", $extParaName[$k]);
	    $tval = $extParaArray[($count%$nrep)*$nName+$k];
	    print FIN "$tname$tval\n";
	}
	print FIN "dxmax               $dxmax\n";
	print FIN "dhmax               $dhmax\n\n";
	print FIN "rStep               $rStep\n";
	print FIN "eStep               $eStep\n";
	print FIN "hStep               $hStep\n";
	print FIN "showHeader          1\n";
	print FIN "dAction             $dAction\n";
	print FIN "dStateStep          $dStateStep\n";
	print FIN "dStatStep           $dStatStep\n";
	print FIN "dInfoStep           $dInfoStep\n\n";
	print FIN "workdir             $cdir\n";
	print FIN "actionFile          $actionFile\n";
	print FIN "stateFile           $stateFile\n";
	print FIN "statFile            $statFile\n";
	print FIN "infoFile            $infoFile\n";
	print FIN "digits              $digits\n\n";
	close(FIN) || die "ERROR: Cannot close $iFile: $!\n";
	system("cat $mFile >> $iFile");
	system("echo '' >> $iFile");
	print "./run$simTag < $iFile > $jFile\n\n";
    }
}

# repair mode: repair possibly bad MC log files
if($task =~ m/^repair\s+([0-9]+)\s+([0-9]+)/i) {
    $nmols = $1;
    $freq = $2;
    $cmd .= "";
    foreach $cdir (@dirArray) {
	$iFile = ($cdir)."/".($inFile);
	$iFileTemp = ($iFile).($tag);
	$cmd .= "./repairMCLogs.pl ";
	$cmd .= ($cdir)."/".($inFile)." ";
	$cmd .= ($cdir)."/".($actionFile)." ";
	$cmd .= ($cdir)."/".($stateFile)." ";
	$cmd .= ($cdir)."/".($statFile)." ";
	$cmd .= ($cdir)."/".($infoFile)." ";
	$cmd .= ($nmols)." ";
	$cmd .= ($freq)." ";
	$cmd .= ($digits)." ";
	$cmd .= ($tag)."\n\n";
    }
    print "$cmd";
}

# resume mode: 
if($task =~ m/^resume\s+([0-9]+)\s+([0-9]+)/i) {
    $rStep = $1;
    $eStep = $2;
    $xhflag = 1;
    $cmd .= "";
    foreach $cdir (@dirArray) {
	$seed = int(rand(1.0e9));
	$cmd .= "./resumeMC.pl ";
	$cmd .= ($seed)." ";
	$cmd .= ($rStep)." ";
	$cmd .= ($eStep)." ";
	$cmd .= ($cdir)."/".($inFile)." ";
	$cmd .= ($cdir)."/".($stateFile)." ";
	$cmd .= ($cdir)."/".($infoFile)." ";
	$cmd .= ($xhflag)." ";
	$cmd .= ($tag)."\n\n";
    }
    print "$cmd";
}

# anal mode: analyze stat.log and state.log
if($task =~ m/^(anal.*)\s+([0-9]+)\s+([0-9]+)/i) {
    $analMols = $1;
    $skipMols = $2;
    $drecMols = $3;
    $beta = $dStateStep/$dStatStep;
    $cmd .= "";
    foreach $cdir (@dirArray) {
	$cmd .= "./anal$simTag ";
	$cmd .= ($cdir)."/".($statFile)." ";
	$cmd .= ($skipMols*$beta)." ";
	$cmd .= ($drecMols*$beta)." ";
	$cmd .= ($digits)." ";
	$cmd .= "> ";
	$cmd .= ($cdir)."/anal.".($statFile)."; ";
	$cmd .= "./$analMols ";
	$cmd .= ($cdir)."/".($stateFile)." ";
	$cmd .= ($skipMols)." ";
	$cmd .= ($drecMols)." ";
	$cmd .= ($digits)." ";
	$cmd .= "> ";
	$cmd .= ($cdir)."/anal.".($stateFile)."\n\n";
    }
    print "$cmd";
}
