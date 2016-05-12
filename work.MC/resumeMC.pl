#! /usr/bin/perl

$nargs = $#ARGV+1;
if($nargs != 7 && $nargs != 8) {
    print "Usage: resumeMC.pl <seed> <rStep> <eStep> <in.log> <state.log> <info.log> <dxmax||dxmax+dhmax> [tag]\n";
    exit(1);
}
$seed = $ARGV[0];
$rStep = $ARGV[1];
$eStep = $ARGV[2];
$inFile = $ARGV[3];
$stateFile = $ARGV[4];
$infoFile = $ARGV[5];
$xhflag = $ARGV[6];
if($nargs == 8) {
    $tag = $ARGV[7];
}
else {
    $tag = ".temp";
}

## main

$inFileTemp = ($inFile).($tag);
$stateFileTemp = ($stateFile).($tag);

$tline = `tail -n3 $infoFile`;
@tinfo = split(/\s+/,$tline);
if($xhflag) {
    $infoColumns = 7;
}
else {
    $infoColumns = 5;
}

open(FTEMP, ">$inFileTemp") || die "ERROR: Cannot open $inFileTemp: $!\n";
open(FIN,   "<$inFile")   || die "ERROR: Cannot open $inFile: $!\n";
while($tline = <FIN>) {
    if($tline =~ m/^(seed)(\s+)/i) {
	if($seed >= 0) {
	    print FTEMP "$1$2$seed\n";
	}
	else {
	    print FTEMP $tline;   
	}
    }
    elsif($tline =~ m/^(dxmax)(\s+)/i) {
	if($#tinfo+1 >= $infoColumns) {
	    $dxmax = $tinfo[$#tinfo-2];
	    print FTEMP "$1$2$dxmax\n";
	}
	else {
	    print FTEMP $tline;
	}
    }
    elsif($xhflag && $tline =~ m/^(dhmax)(\s+)/i) {
	if($#tinfo+1 >= $infoColumns) {
	    $dhmax = $tinfo[$#tinfo-1];
	    print FTEMP "$1$2$dhmax\n";
	}
	else {
	    print FTEMP $tline;
	}
    }
    elsif($tline =~ m/^(rStep)(\s+)/i) {
	print FTEMP "$1$2$rStep\n";
    }
    elsif($tline =~ m/^(eStep)(\s+)/i) {
	print FTEMP "$1$2$eStep\n";
    }
    elsif($tline =~ m/^(showHeader)(\s+)/i) {
	print FTEMP "$1$2"."0\n";
    }
    elsif($tline =~ m/^Energy/) {
	last;
    }
    else {
	print FTEMP $tline;
    }
}
close(FIN)   || die "ERROR: Cannot close $inFile: $!\n";
close(FTEMP) || die "ERROR: Cannot close $inFileTemp: $!\n";

`wc -l $inFile` =~ m/\s*([0-9]+)\s*/i;
$num = $1;
system("tail -n ".(2*$num)." $stateFile > $stateFileTemp");
$result = `./extract -ka $stateFileTemp [-1] -o $inFileTemp`;
system("echo '' >> $inFileTemp");

system("cat $inFile >> ".($inFile).".history");
system("mv -f $inFileTemp $inFile");
system("rm -f $stateFileTemp");

