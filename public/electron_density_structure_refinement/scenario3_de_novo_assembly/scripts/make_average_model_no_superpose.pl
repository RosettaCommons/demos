#!/usr/bin/perl
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
#  @author Frank DiMaio
#  @author Ray Yu-Ruei Wang, wangyr@u.washington.edu
#
##
##
###############################################################################

use strict;
use List::Util qw[min max];

if ($#ARGV < 0) {
	print STDERR "usage: $0 <pdb>*\n";
	exit -1;
}

my %models;
my %lines;
my @idlist;
my ($minres,$maxres) = (9999,0);

foreach my $pdbfile (@ARGV) {
	open (PDB, $pdbfile) || print STDERR "Cannot open $_";
	while (<PDB>) {
		if (! /^ATOM/) {
			next;
		}
		my $atom = substr ($_, 13, 3);
		my $chain = substr($_, 21, 1);
		my $res = int( substr($_, 22, 4) );
		my $id = $res.$atom;
		$minres = min( $minres, $res );
		$maxres = max( $maxres, $res );
	
		my $x = substr ($_, 30, 8);
		my $y = substr ($_, 38, 8);
		my $z = substr ($_, 46, 8);
	
		push @{ $models{ $id } }, [$x,$y,$z];
		$lines{ $id } = $_;
	}
	close (PDB);
}

my $outmean = "average.pdb";

open (OUT1, ">$outmean") || print STDERR "Cannot open $_";
foreach my $resid ( $minres .. $maxres ) {
	foreach my $atomid ("N  ","CA ","C  ","O  ") {
		next if (!defined $models{ $resid.$atomid });

		my $sumX=[0,0,0];
		my $sumX2=[0,0,0];
		my $N = scalar @{ $models{ $resid.$atomid } };
		next if ($N<1);
	
		foreach my $coord (@{ $models{ $resid.$atomid } }) {
			foreach my $i (0..2) {
				$sumX->[$i] += $coord->[$i];
				$sumX2->[$i] += $coord->[$i]*$coord->[$i];
			}
		}
		foreach my $i (0..2) {
			$sumX->[$i] /= $N;
			$sumX2->[$i] = $sumX2->[$i]/$N - $sumX->[$i]*$sumX->[$i];
		}
		#my $B = 0.01 * 8/3 * (3.1415) * (3.1415) * ($sumX2->[0] + $sumX2->[1] + $sumX2->[2]);
        #my $B = sqrt($sumX2->[0] + $sumX2->[1] + $sumX2->[2]);
		#if ($B>99.99) { $B=99.99; }
	
		my $average_line = $lines{ $resid.$atomid };
	
		substr ($average_line, 30, 8) = sprintf ("%8.3f", $sumX->[0]);
		substr ($average_line, 38, 8) = sprintf ("%8.3f", $sumX->[1]);
		substr ($average_line, 46, 8) = sprintf ("%8.3f", $sumX->[2]);
        #substr ($average_line, 60, 7) = sprintf ("%7.3f", $B);
	
		print OUT1 $average_line;
	}
}

close (OUT1);
close (OUT2);




exit 0;
