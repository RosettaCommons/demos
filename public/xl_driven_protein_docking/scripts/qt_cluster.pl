#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 28.01.2013

use strict;
use warnings;
use Getopt::Long;

##############################################################################
##############################################################################
### Extracts the largest clusters from a Quality Threshold matrix.         ###
##############################################################################
##############################################################################
### This script is free of charge for academical use. It is  provided      ###
### WITHOUT warranty of any kind.                                          ###
##############################################################################

my(
    $help,
    $matrix,
    $pdbDir,
    $verbose,
);

&GetOptions(
   "help!"    => \$help,    # prints this help.
   "verbose!" => \$verbose, # outputs additional information on STDERR.
   "infile=s" => \$matrix,  # path to an all against all QT matrix.
   "dir=s"    => \$pdbDir,  # path to the directory that holds the PDB files from the matrix.
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

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

    $usage = "cat xxx.pdb | $0 -o B\n";


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

sub qt_cluster{
    (my $qtMatrix, my $pdbDir) = @_;
 
    my $output = "";

    open(F, $qtMatrix) or die "Failed to open $qtMatrix $!\n";
    my %all;
    my $i = 0; 
    while(<F>){
	chomp($_);
	$all{$i} = $_;
	$i++; 
    }
    close(F);

    my $clusterCenterCount = 0;
    my %remain = %all;
    while((keys %remain) > 0){
	# find decoy with most similar decoys.
	my $max = -1;
	my $maxI = -1;
	my $maxId = "";
	print STDERR "Considering ".(keys %remain)." entries for clustering.\n"
	    if(defined $verbose);

	foreach $i (keys %remain){
	    (my $id, my @a) = split(/\s/, $remain{$i});
	    next unless(-e $id);
	    my $n=0;
	    foreach my $a (@a){
		$n++ if($a==0);
	    }
	    print STDERR "$id:$n vs $maxId:$max\n" if(defined $verbose);
	    if($n > $max){
		$max = $n;
		$maxI = $i;
		$maxId = $id;
	    }
	    elsif($n==$max){
		my $maxEnergy = &getPoseScore("$pdbDir/$maxId");
		my $curEnergy = &getPoseScore("$pdbDir/$id");

		chomp($maxEnergy);
		chomp($curEnergy);

		if($curEnergy < $maxEnergy){
		    print STDERR "$id and maxId $maxId have both $n number of cluster members. ".
			"Prefering $id with lower score $curEnergy ($maxId: $maxEnergy).\n"
			if(defined $verbose);
		    $max = $n;
		    $maxI = $i;
		    $maxId = $id;
		}
	    }
	}
	print STDERR "Using $maxI. $maxId as cluster center\n" if(defined $verbose);
	$clusterCenterCount++;
	$output .= "$maxId\t$clusterCenterCount\n";
	print STDERR "remove decoy and all similar decoys from remaining list by row\n"
	    if(defined $verbose);

	my $col2del = "#$maxI#";
	(my $id, my @a) = split(/\s/, $remain{$maxI});
	delete($remain{$maxI});

	for($i=0; $i<=$#a; $i++){
	    if($a[$i] == 0){
		unless (exists $remain{$i}){
		    print STDERR "Have removed cluster center itself: $i $id\n" if(defined $verbose);
		    next;
		}
		print STDERR "Removing $i. ".substr($remain{$i},0,index($remain{$i}," ")).
		      " from cluster list, as it is a cluster member to $maxId\n"
		      if(defined $verbose);
		delete($remain{$i});
		$col2del .= "#$i#";
	    }
	}

	print STDERR "Now delete all entries that are not similar to any other entry, which will prohibit ".
	      "them from being assigned to a cluster center.\n"
	      if(defined $verbose);
	foreach my $j (sort keys %remain){
	    (my $id, my @b) = split(/\s/, $remain{$j});
	    my $hasSim = 0;
	    for($i=0; $i<=$#b; $i++){
		if($b[$i] == 0 and $i != $j){
		    $hasSim = 1;
		}
	    }
#	    if($hasSim == 0){
#		print STDERR "Removing $j. $id from cluster list, as it has no similarity to remaining clusters.\n".
#		    $remain{$j}."\n" if(defined $verbose);
#		delete($remain{$j});
#		$col2del .= "#$j#";
#	    }
	}
	print STDERR "Following columns will be deleted: $col2del\n" if(defined $verbose);

	print STDERR "Now deleting columns in allVsAll RMSD matrix of just deleted entries\n".
	    "Remaining ".(keys %remain)." entries after line-wise cluster member removement.\n"
	    if(defined $verbose);
	foreach my $x (sort keys %remain){
	    print STDERR "$remain{$x}\n" if(defined $verbose);
	}
	foreach my $j (sort keys %remain){
	    (my $id, my @b) = split(/\s/, $remain{$j});
	    # no need to explicitly delete column of maxI as it is automatically
	    # included in @a due to self comparison = 0;
	    my $newLine = "";
	    for($i=0; $i<=$#b; $i++){
		# deletion done by simply not adding RMSD value back to remaining matrix
		next if($col2del =~ /\#$i\#/);
		# undeleted line should be dissimilar to the deleted entries.
		print STDERR "WARNING $i!\n" if($a[$i] == 0 and defined $verbose);
		$newLine .= "$b[$i] ";
	    }
	    # remove last space character
	    chop($newLine);
	    if(length($newLine)!=0){
		$remain{$j} = "$id $newLine";
	    }
	}

	my %remain2=();
	$i=0;
	foreach my $x (sort {$a <=> $b} keys %remain){
	    $remain2{$i} = $remain{$x}; 
	    $i++
	}
	%remain = %remain2;

	print STDERR "Remaining ".(keys %remain)." entries after line-wise + column-wise cluster member deletion.\n"
	    if(defined $verbose);
	foreach my $x (sort { $a<=>$b } keys %remain){
	    print STDERR "$remain{$x}\n" if(defined $verbose);
	}
    }

    return $output;
}

 ##############################################################################

sub getPoseScore {
    (my $fileName) = @_;
    my $command = "grep I_sc $fileName | sed 's/.* //'"; 
    print STDERR "Extracting interface score: $command.\n" if(defined $verbose);
    my $e = `$command`;
    if(length($e)==0){
        $command = "grep ^pose $fileName | sed 's/.* //'"; 
	print STDERR "Interface score not found. Extracting pose score: $command.\n"
	    if(defined $verbose);
	$e = `$command`;
    }
    chomp($e);
    print STDERR "Score is: $e\n" if(defined $verbose);
    return $e;
}

##############################################################################
##############################################################################
### MAIN
##############################################################################
##############################################################################

if(!defined $matrix or !defined $pdbDir){
    die "\nPlease specify all required options. Try \"$0 -h\" for a complete list of options\n\n";
} else {
    die "\nCan't find the file $matrix. Please ensure that it exists!\n\n" if(!-e $matrix);
    die "\nCan't find the directory $pdbDir. Please ensure that it exists!\n\n" if(!-d $pdbDir);
    
    my $cluster = &qt_cluster($matrix, $pdbDir);
    print $cluster;
}

