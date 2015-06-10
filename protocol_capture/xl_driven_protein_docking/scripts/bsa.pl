#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 23.01.2013

use strict;
use warnings;
use Cwd;
use File::Basename;
use Getopt::Long;
use File::Path;

##############################################################################
### Calculates the buried solvent area between all components of a protein ###
### complex.                                                               ###
##############################################################################
##############################################################################
### This script is free of charge for academical use. It is  provided      ###
### WITHOUT warranty of any kind.                                          ###
##############################################################################

my $tmpDir = "Trash";
my $tmpOut = "2bDeleted";

my (
    # variable for parameters which are read in from commandline
    $help,
    $naccessPath,
    $verbose,
    $xDelTrash,
   );

##############################################################################
### read all needed parameters from commandline ##############################
##############################################################################

&GetOptions(
    "help!"   => \$help,        # prints this help.
    "naccess=s" => \$naccessPath, # path to the NACCESS executable
    "v!"      => \$verbose,     # prints out additional information
    "d!"      => \$xDelTrash,   # does not delete temporary directory Trash
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################
# help
if ($help) {printHelp(); exit}

##############################################################################
##############################################################################
### SUBROUTINES
##############################################################################
##############################################################################
###############################################################################
sub printHelp {
###############################################################################
    # prints a help about the using and parameters of this scripts 
    # (execute if user types commandline parameter -h)
    # param:  no paramaters
    # return: no return value

    my (
	$usage,
	$sourceCode,
	@rows,
	$row,
	$option,
	$scriptInfo,
	$example,
       );

    $usage = "$0 xxx.pdb -naccess /path/to/naccess/executable\n";


    print "\nUsage: " .  $usage . "\n";

    print "Valid options are:\n\n";
    open(MYSELF, "$0") or
      die "Cannot read source code file $0: $!\n";
    $sourceCode .= join "", <MYSELF>;
    close MYSELF;
    $sourceCode =~ s/^.+?\&GetOptions\(\n//s;
    $sourceCode =~ s/\n\).+$//s;
    @rows = split /\n/, $sourceCode;
    foreach $row (@rows){
        $option = $row;
	$option =~ s/\s+\"//g;
	$option =~ s/\"\s.+\#/\t\#/g;
	$option =~ s/=./\t<value> [required]/;
	$option =~ s/:./\t<value> [optional]/;
	$option =~ s/!/\t<non value> [optional]/;

	$row =~ s/^.*//;
	print "\t";
	printf("%-1s%-30s%-30s\n", "-",$option,$row);

    } # end of foreach $row (@rows)
    print "\n";
    print "Options may be abreviated, e.g. -h for --help\n\n";

    $example  = "$0";
}
##############################################################################
# run NACCESS and get total SAS on PDB file.
sub naccess{

    (my $fileName) = @_;

    my $command = "$naccessPath $fileName.pdb"; 
    if(defined $verbose){
	print STDERR "Executing $command\n";
	system($command);
    }
    else{
	`$command`;
    }
    if(-e "$fileName.asa"){
	$command = "grep TOTAL $fileName.rsa";
	print STDERR "Executing $command\n" if(defined $verbose);
	(my $flag, my $totalSas, my $sideChainSas, my $mainChainSas, my $polarSas, my $nonPolarSas)
	    = split(/\s+/,`$command`);
 
	return $totalSas;
    }
}
##############################################################################
# returns all chains from a PDB file.
sub getNumberOfPDBchains {
    open(F,$ARGV[0]) or die "Couldn't open file \"$ARGV[0]\"\n";
    my %h;
    while(<F>){
	if(/^ATOM/){
	    my $chainId  = substr($_, 21, 1);
	    $chainId = "_" if($chainId eq " ");
	    $h{$chainId}="";
	}
    }
    close(F);
    my @chains = ();
    foreach my $c (sort keys %h){
	push(@chains, $c);
    }
    return \@chains;
}
##############################################################################
# extracts certain chains from a PDB file.
sub extractChainIntoFile {

    (my $infile, my $chainId, my $outfile) = @_;

    open(F,$infile) or die "Couldn't open file \"$infile\"\n";
    open(O,">$outfile") or die "Couldn't open file \"$outfile\"\n";

    while(<F>){
	if(/^ATOM/){
	    my $c = substr($_, 21, 1);
	    if($c eq $chainId){
		print O $_;
	    }
	}
    }
    close(F);
    close(O);
}

##############################################################################
##############################################################################
### MAIN
##############################################################################
##############################################################################

if (defined $ARGV[0] and defined $naccessPath) {
    if(-e $ARGV[0]){
	my $complex = getcwd()."/".$ARGV[0];
	my @chains = @{&getNumberOfPDBchains()};

	print STDERR "Creating directory $tmpDir\n" if(defined $verbose);
	mkdir($tmpDir);
	print STDERR "Changing to directory $tmpDir\n" if(defined $verbose);
	chdir($tmpDir);

	for(my $i=0; $i<=$#chains; $i++){
	    for(my $j=$i+1; $j<=$#chains; $j++){

		my $c1 = $chains[$i];
		&extractChainIntoFile($complex, $c1, "$tmpOut-$c1.pdb");
		my $sas1 = &naccess("$tmpOut-$c1");
		
		my $c2 = $chains[$j];
		&extractChainIntoFile($complex, $c2, "$tmpOut-$c2.pdb");
		my $sas2 = &naccess("$tmpOut-$c2");

		system("cat $tmpOut-$c1.pdb $tmpOut-$c2.pdb > $tmpOut-$c1$c2.pdb");
		my $sas12 = &naccess("$tmpOut-$c1$c2");

		my $sbsa = $sas1 + $sas2 - $sas12;
		print STDERR "file\tchainId1\tchainId2\tsbsa\tsas1\tsas2\tsas1\&2\n" if(defined $verbose);
		printf("%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\n",basename($ARGV[0]),$c1,$c2,$sbsa,$sas1,$sas2,$sas12);
	    }
	}

	print STDERR "Changing back to starting directory\n" if(defined $verbose);
	chdir("../");
	if(!defined $xDelTrash){
	    print STDERR "Removing temporary directory $tmpDir\n" if(defined $verbose);
	    rmtree($tmpDir);
	}
    }
}
else{
    die "\nUsage: $0 complex.pdb\n\n";
}
