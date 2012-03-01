#!/usr/bin/perl
##
##
## Copyright 2003, University of Washington, the Baker Lab, and Dylan Chivian.
##   This document contains private and confidential information and
##   its disclosure does not constitute publication.  All rights are
##   reserved by University of Washington, the Baker Lab, and Dylan Chivian,
##   except those specifically granted by license.
##
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.1.1.1 $
##  $Date: 2003/09/05 01:47:28 $
##  $Author: dylan $
##
##
###############################################################################

###############################################################################
package PDButil;
###############################################################################

###############################################################################
# conf
###############################################################################

$| = 1;                                                   # don't buffer stdout

$cst_radius = 8.0;  # in angstroms
$cst_pad    = 0.0;  # in angstroms

###############################################################################
# init
###############################################################################

local %opts = &getCommandLineOptions ();
local $fastafile           = $opts{fastafile};
local $parent_pdb          = $opts{parentpdb};
local $parent_ssa          = $opts{parentssa};
local $zonesfile           = $opts{zonesfile};
local $outpdb              = $opts{outpdb};
local $outssa              = $opts{outssa};
local $outcst              = $opts{outcst};
#local $parent_fssp        = $opts{parentfssp};
local $takeoffpad_flag     = $opts{takeoffpad};
local $sidechains_flag     = $opts{sidechains};
local $loopregiononly_flag = $opts{loopregiononly};
local $keep_hetero         = $opts{keephetero};

###############################################################################
# main
###############################################################################

# read query fasta
#
@query_fasta = ();
@query_fasta_buf = &fileBufArray ($fastafile);
foreach $line (@query_fasta_buf) {
    next if ($line =~ /^\>/);
    $line =~ s/\s+//g;
    push (@query_fasta, split (//, $line));
}
$seq_len = $#query_fasta + 1;


# read parent pdb
#
@pdb_buf = &fileBufArray ($parent_pdb);
foreach $line (@pdb_buf) {
    if ($line =~ /^ATOM/) {
	$atomtype          = substr ($line, 12, 4);
	$restype           = substr ($line, 17, 3);
	$res_i             = substr ($line, 22, 4) - 1;

	# store parent sequence so we can check for identical residues
	$parent_fasta[$res_i] = &mapResCode ($restype);

	# need occupancy of parent, use CA occ
	$p_occ[$res_i] = substr ($line, 54, 6)  if ($atomtype eq ' CA ');

	# N, CA, C, O, CB only
	if (defined &typeI($atomtype)) {
	    $src->[$res_i]->[&typeI($atomtype)]->[0] = substr ($line, 30, 8);
	    $src->[$res_i]->[&typeI($atomtype)]->[1] = substr ($line, 38, 8);
	    $src->[$res_i]->[&typeI($atomtype)]->[2] = substr ($line, 46, 8);
	}

	# handle template from glycine (don't worry: will be replaced if CB available)
	if (!  defined $src->[$res_i]->[&typeI('CB')]
	    && defined $src->[$res_i]->[&typeI('N')]
	    && defined $src->[$res_i]->[&typeI('CA')]
	    && defined $src->[$res_i]->[&typeI('C')]
	    && $p_occ[$res_i] > 0
	   ) {
	    $src->[$res_i]->[&typeI('CB')] = 
		&getCbCoords ($src->[$res_i]->[&typeI('N')],
			      $src->[$res_i]->[&typeI('CA')],
			      $src->[$res_i]->[&typeI('C')],
			      $restype);
	}
	
	# store parent sequence and all atoms so we can recover side-chains
	$all_atoms->[$res_i]->{$atomtype}->[0] = substr ($line, 30, 8);
	$all_atoms->[$res_i]->{$atomtype}->[1] = substr ($line, 38, 8);
	$all_atoms->[$res_i]->{$atomtype}->[2] = substr ($line, 46, 8);
    }
    elsif ($line =~ /^HETATM/ && $line !~ /HOH/) {
	push (@hetatms, $line);
    }
}
# debug
#for ($pj=0; $pj <= $#p_occ; ++$pj) {
#    print "p_occ[$pj] = $p_occ[$pj]\n";
#}


# read zones file
#
@zones_buf = &fileBufArray ($zonesfile);
foreach $line (@zones_buf){
    if ($line =~ /^\s*zone\s*\:?\s*none\s*/i) {
	$empty_correct = 'true';
	last;
    }
    next if ($line !~ /^\s*zone\s*\:?\s*(\d+)\s*\-\s*(\d+)\s*\:\s*(\d+)\s*\-\s*(\d+)\s*/i);
    $q_start = $1;
    $q_stop  = $2;
    $p_start = $3;
    $p_stop  = $4;
    if ($loopregiononly_flag !~ /^t/i) {
	if ($q_stop - $q_start != $p_stop - $p_start) {
	    &abort ("unequal zone $q_start-$q_stop:$p_start-$p_stop");
	}
    }
    --$q_start;  --$q_stop;  --$p_start;  --$p_stop;

# debug
#print "$q_start-$q_stop:$p_start-$p_stop\n";

    # for loop regions, only care about edges of zones
    #
    if ($loopregiononly_flag =~ /^t/i) {
	# stem for Cterm loop
	if ($q_start == 0) {
	    # last aligned
	    &abort ("query is aligned to missing density for query[$q_stop] with parent[$p_stop]")  if ($p_occ[$p_stop] <= 0);
	    $q2p_mapping[$q_stop] = $p_stop;

#	    # aligned adjacent to last aligned
#	    if ($q_stop - $q_start > 0) {
#		$q_adj = $q_stop-1;
#		$p_adj = $p_stop-1;
#		&abort ("query is aligned to missing density for query[$q_adj] with parent[$p_adj]")  if ($p_occ[$p_adj] <= 0);
#		$q2p_mapping[$q_adj] = $p_adj;
#	    }
	}
	# stem for Nterm loop
	elsif ($q_stop == $seq_len - 1) {
	    # last aligned
	    &abort ("query is aligned to missing density for query[$q_start] with parent[$p_start]")  if ($p_occ[$p_start] <= 0);
	    $q2p_mapping[$q_start] = $p_start;

#	    # aligned adjacent to last aligned
#	    if ($q_stop - $q_start > 0) {
#		$q_adj = $q_start+1;
#		$p_adj = $p_start+1;
#		&abort ("query is aligned to missing density for query[$q_adj] with parent[$p_adj]")  if ($p_occ[$p_adj] <= 0);
#		$q2p_mapping[$q_adj] = $p_adj;
#	    }
	}
	# stems for internal loop
	else {
	    # last aligned
	    &abort ("query is aligned to missing density for query[$q_start] with parent[$p_start]")  if ($p_occ[$p_start] <= 0);
	    &abort ("query is aligned to missing density for query[$q_stop] with parent[$p_stop]")  if ($p_occ[$p_stop] <= 0);
	    $q2p_mapping[$q_start] = $p_start;
	    $q2p_mapping[$q_stop] = $p_stop;

#	    # aligned adjacent to last aligned
#	    if ($q_stop - $q_start > 0) {
#		# start adj
#		$q_adj = $q_start+1;
#		$p_adj = $p_start+1;
#		&abort ("query is aligned to missing density for query[$q_adj] with parent[$p_adj]")  if ($p_occ[$p_adj] <= 0);
#		$q2p_mapping[$q_adj] = $p_adj;
#
#		# stop adj
#		$q_adj = $q_stop-1;
#		$p_adj = $p_stop-1;
#		&abort ("query is aligned to missing density for query[$q_adj] with parent[$p_adj]")  if ($p_occ[$p_adj] <= 0);
#		$q2p_mapping[$q_adj] = $p_adj;
#	    }
	}
    }

    # for full template, do full body of zone
    #
    else {
	if ($p_start >= 9998 || $p_stop >= 9998) {
	    &abort ("no available parent coords");
	}
	for ($i=0; $i <= $q_stop-$q_start; ++$i) {
	    $qi = $q_start + $i;
	    $pj = $p_start + $i;
	    &abort ("query is aligned to missing density for query[$qi] with parent[$pj]")  if ($p_occ[$pj] <= 0);

	    $q2p_mapping[$qi] = $pj;
	}
    }
}
if ($#q2p_mapping < 0) {
    if ($empty_correct) {
	print STDERR "zone NONE in $zonesfile\n";
	exit 0;
    }
    &abort ("no zones read in $zonesfile\n");
}
# debug
#for ($qi=0; $qi <= $#q2p_mapping; ++$qi) {
#    $pj = $q2p_mapping[$qi];
#    print "MAP '$qi' '$pj'\n";
#}



# build template ssa and output
#
if ($outssa) {

    # read parent ssa
    #
    @parent_ssa = ();
    @parent_ssa_buf = &fileBufArray ($parent_ssa);
    foreach $line (@parent_ssa_buf) {
	next if ($line =~ /^\>/);
	$line =~ s/\s+//g;
	push (@parent_ssa, split (//, $line));
    }

    # assign query ssa (L for unaligned regions)
    #
    $outssa_str = '';
    for ($qi=0; $qi <= $#query_fasta; ++$qi) {
	$pj = $q2p_mapping[$qi];

	if (defined $pj) {
	    $query_ssa[$qi] = ($parent_ssa[$pj]) ? $parent_ssa[$pj] : 'L';
	}
	else {
	    $query_ssa[$qi] = 'L';
	}
	$outssa_str .= $query_ssa[$qi];
	#$outssa_str .= "\n"  if ($qi % 50 == 0);
    }

    # supress short regions of regular ss which are likely wrong
    $outssa_str =~ s/^(H{1,2})([^H])/('L'x(length($1))).$2/ei;
    $outssa_str =~ s/([^H])(H{1,2})$/$1.('L'x(length($2)))/ei;
    $outssa_str =~ s/([^H])(H{1,2})([^H])/$1.('L'x(length($2))).$3/eig;
    $outssa_str =~ s/^(E{1})([^E])/('L'x(length($1))).$2/ei;
    $outssa_str =~ s/([^E])(E{1})$/$1.('L'x(length($2)))/ei;
    $outssa_str =~ s/([^E])(E{1})([^E])/$1.('L'x(length($2))).$3/eig;

    if ($outssa) {
	open (OUTSSA, '>'.$outssa);
	print OUTSSA "$outssa_str\n";
	close (OUTSSA);
    }
}


# build distance constraints
#
if ($outcst) {

    # measure parent pdb distances
    #
    for ($p_i=0; $p_i <= $#{$src}; ++$p_i) {
	$p_dist->[$p_i]->[$p_i]->{CB}->{CB} = 0;
	for ($p_j=$p_i+1; $p_j <= $#{$src}; ++$p_j) {
	    $p_dist->[$p_i]->[$p_j]->{CB}->{CB}
	    = $p_dist->[$p_j]->[$p_i]->{CB}->{CB} = &measureDist ($src->[$p_i], 
								  $src->[$p_j],
								  'CB',
								  'CB');
	}
    }

    # map to query distances
    #
    @constraints = ();
    for ($q_i=0; $q_i <= $#query_fasta; ++$q_i) {
	for ($q_j=$q_i+5; $q_j <= $#query_fasta; ++$q_j) {
	    
	    if (defined $q2p_mapping[$q_i] && defined $q2p_mapping[$q_j]) {
		$p_i = $q2p_mapping[$q_i];
		$p_j = $q2p_mapping[$q_j];
		
		if ($p_occ[$p_i] > 0 && $p_occ[$p_j] > 0) {
		    $dist = $p_dist->[$p_i]->[$p_j]->{CB}->{CB};
		    if ($dist <= $cst_radius) {
			$dist += $cst_pad;
			push (@constraints, 
			      sprintf (" %6d %2s %6d %2s   %10.2f %10.2f %10.2f",
				       $q_i+1, 'CB', $q_j+1, 'CB', $dist, 0.00, 0.00));
		    }
		}
	    }
	}
    }    

    # write cst file
    #
    $num_recs = $#constraints + 1;
    @constraints = ("NMR_v3.0", 
		    "CB-CB csts from ksync alignment", 
		    "associated with $outpdb",
		    $num_recs,
		    @constraints
		    );
    open  (OUTCST, '>'.$outcst);
    print  OUTCST join ("\n", @constraints) . "\n";
    close (OUTCST);
}


# build template pdb
#
@atoms = ();
for ($qi=0; $qi <= $#query_fasta; ++$qi) {
    $q_restype = $query_fasta[$qi];
    if ($sidechains_flag =~ /^F/i) {
	@atom_recs = &resAtomRecsNoSidechain ($q_restype);
    } else {
	@atom_recs = &resAtomRecs ($q_restype);
    }
    for ($i=0; $i <= $#atom_recs; ++$i) {
	substr ($atom_recs[$i], 22, 4) = sprintf ("%4d", $qi+1);      # res_num
    }
    if (defined $q2p_mapping[$qi] ||
	($takeoffpad_flag =~ /^T/i && ! defined $q2p_mapping[$qi] && 
	 (($qi != 0 && defined $q2p_mapping[$qi-1]) || defined $q2p_mapping[$qi+1])
	)
       ) {

	if (defined $q2p_mapping[$qi]) {
	    $pj = $q2p_mapping[$qi];
	} elsif ($qi != 0 && defined $q2p_mapping[$qi-1]) {
	    $pj = $q2p_mapping[$qi-1]+1;
	} else {
	    $pj = $q2p_mapping[$qi+1]-1;
	}

	if (defined $q2p_mapping[$qi] && $p_occ[$pj] <= 0) {
	    &abort ("query aligned to missing density from query[$qi] to parent[$pj]");
	} elsif ($p_occ[$pj] <= 0) {
	    &abort ("attempt to assign missing density from parent[$pj] as a takeoff for query[$qi]");
	}

	# replace CB if not identical residue (distances vary)
	if ($query_fasta[$qi] ne $parent_fasta[$pj]) {
	    $src->[$pj]->[&typeI('CB')] =
		&getCbCoords ($src->[$pj]->[&typeI('N')],
			      $src->[$pj]->[&typeI('CA')],
			      $src->[$pj]->[&typeI('C')],
			      $q_restype);
	} 

	# fill backbone coords
	for ($i=0; $i < 5 && defined $atom_recs[$i]; ++$i) {
	    substr ($atom_recs[$i], 54, 6) = sprintf ("%6.2f", 1.00);
	    substr ($atom_recs[$i], 30, 8) = sprintf ("%8.3f", $src->[$pj]->[$i]->[0]);
	    substr ($atom_recs[$i], 38, 8) = sprintf ("%8.3f", $src->[$pj]->[$i]->[1]);
	    substr ($atom_recs[$i], 46, 8) = sprintf ("%8.3f", $src->[$pj]->[$i]->[2]);
	}

	# add side-chains if identical residue
	if ($sidechains_flag !~ /^F/i) {
	    if ($query_fasta[$qi] eq $parent_fasta[$pj]) {
		for ($i=5; $i <= $#atom_recs; ++$i) {
		    $atomtype = substr ($atom_recs[$i], 12, 4);
		    if (defined $all_atoms->[$pj]->{$atomtype}->[0]) {
			substr ($atom_recs[$i], 54, 6) = sprintf ("%6.2f", 1.00);
			substr ($atom_recs[$i], 30, 8) = sprintf ("%8.3f", $all_atoms->[$pj]->{$atomtype}->[0]);
			substr ($atom_recs[$i], 38, 8) = sprintf ("%8.3f", $all_atoms->[$pj]->{$atomtype}->[1]);
			substr ($atom_recs[$i], 46, 8) = sprintf ("%8.3f", $all_atoms->[$pj]->{$atomtype}->[2]);
		    }
		}
	    }
	}
    }
    push (@atoms, @atom_recs);
}    


# output
#
@outbuf = (@atoms, "TER    atmi");
if ($keep_hetero) {
    push (@outbuf, @hetatms);
}
for ($i=0; $i <= $#outbuf; ++$i) {
    substr ($outbuf[$i], 6, 5) = sprintf ("%5d", $i+1);
}


# finish template pdb
$outpdb_str = join ("\n", @outbuf) ."\nEND\n";
if ($outpdb) {
#    print "creating $outpdb\n";
    open (OUTPDB, '>'.$outpdb);
    select (OUTPDB);
}
print $outpdb_str;
if ($outpdb) {
    close (OUTPDB);
    select (STDOUT);
}


# exit
#
exit 0;

###############################################################################
# subs
###############################################################################

# getCbCoords()
#
sub getCbCoords {
    my ($N_coords, $Ca_coords, $C_coords, $restype) = @_;
    my $Cb_coords = [];

    # formula (note: all vectors are normalized)
    # CaCb = bondlen * [cos5475*(-CaN -CaC) + sin5475*(CaN x CaC)]

    # config
    my $cos5475 = 0.577145190;                        # cos 54.75 = cos 109.5/2
    my $sin5475 = 0.816641555;                        # sin 54.75 = sin 109.5/2
    my $CC_bond = 1.536;                                          # from ethane
    my %CaCb_bond = ( 'A' => 1.524,
		      'C' => 1.531,
		      'D' => 1.532,
		      'E' => 1.530,
		      'F' => 1.533,
		      'G' => 1.532,
		      'H' => 1.533,
		      'I' => 1.547,
		      'K' => 1.530,
		      'L' => 1.532,
		      'M' => 1.530,
		      'N' => 1.532,
		      'P' => 1.528,
		      'Q' => 1.530,
		      'R' => 1.530,
		      'S' => 1.530,
		      'T' => 1.545,
		      'V' => 1.546,
		      'W' => 1.533,
		      'Y' => 1.533,
                      'ALA' => 1.524,
		      'CYS' => 1.531,
		      'ASP' => 1.532,
		      'GLU' => 1.530,
		      'PHE' => 1.533,
		      'GLY' => 1.532,
		      'HIS' => 1.533,
		      'ILE' => 1.547,
		      'LYS' => 1.530,
		      'LEU' => 1.532,
		      'MET' => 1.530,
		      'ASN' => 1.532,
		      'PRO' => 1.528,
		      'GLN' => 1.530,
		      'ARG' => 1.530,
		      'SER' => 1.530,
		      'THR' => 1.545,
		      'VAL' => 1.546,
		      'TRP' => 1.533,
		      'TYR' => 1.533,
		    );
    my $bondlen = (defined $restype) ? $CaCb_bond{$restype} : $CC_bond;

    # init vectors
    my $CaN  = +[];  my $CaN_mag  = 0.0;
    my $CaC  = +[];  my $CaC_mag  = 0.0;
    my $vert = +[];  my $vert_mag = 0.0;
    my $perp = +[];  my $perp_mag = 0.0;
    my $CaCb = +[];

    # CaN
    for ($i=0; $i<3; ++$i) {
	$CaN->[$i] = $N_coords->[$i] - $Ca_coords->[$i];
	$CaN_mag += $CaN->[$i] * $CaN->[$i];
    }
    $CaN_mag = sqrt ($CaN_mag);
    for ($i=0; $i<3; ++$i) {
	$CaN->[$i] /= $CaN_mag;
    }

    # CaC
    for ($i=0; $i<3; ++$i) {
	$CaC->[$i] = $C_coords->[$i] - $Ca_coords->[$i];
	$CaC_mag += $CaC->[$i] * $CaC->[$i];
    }
    $CaC_mag = sqrt ($CaC_mag);
    for ($i=0; $i<3; ++$i) {
	$CaC->[$i] /= $CaC_mag;
    }

    # vert = -CaN -CaC
    for ($i=0; $i<3; ++$i) {
	$vert->[$i] = - $CaN->[$i] - $CaC->[$i];
	$vert_mag += $vert->[$i] * $vert->[$i];
    }
    $vert_mag = sqrt ($vert_mag);
    for ($i=0; $i<3; ++$i) {
	$vert->[$i] /= $vert_mag;
    }

    # perp = CaN x CaC
    $perp->[0] = $CaN->[1] * $CaC->[2] - $CaN->[2] * $CaC->[1];
    $perp->[1] = $CaN->[2] * $CaC->[0] - $CaN->[0] * $CaC->[2];
    $perp->[2] = $CaN->[0] * $CaC->[1] - $CaN->[1] * $CaC->[0];
    # x product of two unit vectors is already unit, so no need to normalize

    # CaCb
    for ($i=0; $i<3; ++$i) {
	$CaCb->[$i] = $bondlen * ($cos5475 * $vert->[$i] +
				  $sin5475 * $perp->[$i]);
    }

    # Cb_coords
    for ($i=0; $i<3; ++$i) {
	$Cb_coords->[$i] = $Ca_coords->[$i] + $CaCb->[$i];
    }

    return $Cb_coords;
}


# typeI()
#
sub typeI {
    my $atomtype = shift;
    $atomtype =~ s/\s+//g;

    my %typenum = ( 'N'  => 0,
		    'CA' => 1,
		    'C'  => 2,
		    'O'  => 3,
		    'CB' => 4
		  );

    return $typenum{$atomtype};
}


# measureDist()
#
sub measureDist {
    my ($res_1, $res_2, $atomtype_1, $atomtype_2) = @_;
    my $atomtype_1_i = &typeI($atomtype_1);
    my $atomtype_2_i = &typeI($atomtype_2);
    
    my $x_diff = $res_1->[$atomtype_1_i]->[0] - $res_2->[$atomtype_2_i]->[0];
    my $y_diff = $res_1->[$atomtype_1_i]->[1] - $res_2->[$atomtype_2_i]->[1];
    my $z_diff = $res_1->[$atomtype_1_i]->[2] - $res_2->[$atomtype_2_i]->[2];

    return sqrt ($x_diff*$x_diff + $y_diff*$y_diff + $z_diff*$z_diff);
}


# mapResCode()
#
sub mapResCode {
    local ($incode, $silent) = @_;
    $incode = uc $incode;
    my $newcode = undef;

    my %one_to_three = ( 'G' => 'GLY',
			 'A' => 'ALA',
			 'V' => 'VAL',
			 'L' => 'LEU',
			 'I' => 'ILE',
			 'P' => 'PRO',
			 'C' => 'CYS',
			 'M' => 'MET',
			 'H' => 'HIS',
			 'F' => 'PHE',
			 'Y' => 'TYR',
			 'W' => 'TRP',
			 'N' => 'ASN',
			 'Q' => 'GLN',
			 'S' => 'SER',
			 'T' => 'THR',
			 'K' => 'LYS',
			 'R' => 'ARG',
			 'D' => 'ASP',
			 'E' => 'GLU',
			 'X' => 'XXX',
			 '0' => '  A',
			 '1' => '  C',
			 '2' => '  G',
			 '3' => '  T',
			 '4' => '  U'
			);

    my %three_to_one = ( 'GLY' => 'G',
			 'ALA' => 'A',
			 'VAL' => 'V',
			 'LEU' => 'L',
			 'ILE' => 'I',
			 'PRO' => 'P',
			 'CYS' => 'C',
			 'MET' => 'M',
			 'HIS' => 'H',
			 'PHE' => 'F',
			 'TYR' => 'Y',
			 'TRP' => 'W',
			 'ASN' => 'N',
			 'GLN' => 'Q',
			 'SER' => 'S',
			 'THR' => 'T',
			 'LYS' => 'K',
			 'ARG' => 'R',
			 'ASP' => 'D',
			 'GLU' => 'E',
			 
			 '  X' => 'X',
			 '  A' => '0',
			 '  C' => '1',
			 '  G' => '2',
			 '  T' => '3',
			 '  U' => '4',
			 ' +A' => '0',
			 ' +C' => '1',
			 ' +G' => '2',
			 ' +T' => '3',
			 ' +U' => '4',
			 
			 # all of these are supposed to be handled by MODRES
			 #   (but aren't always or aren't done properly)
			 '5HP' => 'Q',
			 'ABA' => 'C',
			 'AGM' => 'R',
			 'CEA' => 'C',
			 'CGU' => 'E',
			 'CME' => 'C',
			 'CSB' => 'C',
			 'CSE' => 'C',
			 'CSD' => 'C',
			 'CSO' => 'C',
			 'CSP' => 'C',
			 'CSS' => 'C',
			 'CSW' => 'C',
			 'CSX' => 'C',
			 'CXM' => 'M',
			 'CYM' => 'C',
			 'CYG' => 'C',
			 'DOH' => 'D',
			 'FME' => 'M',
			 'GL3' => 'G',
			 'HYP' => 'P',
			 'KCX' => 'K',
			 'LLP' => 'K',
			 'LYZ' => 'K',
			 'MEN' => 'N',
			 'MGN' => 'Q',
			 'MHS' => 'H',
			 'MIS' => 'S',
			 'MLY' => 'K',
			 'MSE' => 'M',
			 'NEP' => 'H',
			 'OCS' => 'C',
			 'PCA' => 'Q',
			 'PTR' => 'Y',
			 'SAC' => 'S',
			 'SEP' => 'S',
			 'SMC' => 'C',
			 'STY' => 'Y',
			 'SVA' => 'S',            
			 'TPO' => 'T',
			 'TPQ' => 'Y',
			 'TRN' => 'W',
			 'TRO' => 'W',
			 'YOF' => 'Y',
			 
			 '1MG' => 'X',
			 '2DA' => 'X',
			 '2PP' => 'X',
			 '4SC' => 'X',
			 '4SU' => 'X',
			 '5IU' => 'X',
			 '5MC' => 'X',
			 '5MU' => 'X',
			 'ACB' => 'X',
			 'ACE' => 'X',
			 'ACL' => 'X',
			 'ADD' => 'X',
			 'AHO' => 'X',
			 'AIB' => 'X',
			 'ALS' => 'X',
			 'ARM' => 'X',
			 'ASK' => 'X',
			 'ASX' => 'X',          # NOT B, PREFER TOTAL AMBIGUITY
			 'BAL' => 'X',
			 'BE2' => 'X',
			 'CAB' => 'X',
			 'CBX' => 'X',
			 'CBZ' => 'X',
			 'CCC' => 'X',
			 'CHA' => 'X',
			 'CH2' => 'X',
			 'CH3' => 'X',
			 'CHG' => 'X',
			 'CPN' => 'X',
			 'CRO' => 'X',
			 'DAL' => 'X',
			 'DGL' => 'X',
			 'DOC' => 'X',
			 'DPN' => 'X',
			 'EXC' => 'X',
			 'EYS' => 'X',
			 'FGL' => 'X',
			 'FOR' => 'X',
			 'G7M' => 'X',
			 'GLQ' => 'X',
			 'GLX' => 'X',          # NOT Z, PREFER TOTAL AMBIGUITY
			 'GLZ' => 'X',
			 'GTP' => 'X',
			 'H2U' => 'X',
			 'HAC' => 'X',
			 'HEM' => 'X',
			 'HMF' => 'X',
			 'HPB' => 'X',
			 'IAS' => 'X',
			 'IIL' => 'X',
			 'IPN' => 'X',
			 'LAC' => 'X',
			 'LYT' => 'X',
			 'LYW' => 'X',
			 'MAA' => 'X',
			 'MAI' => 'X',
			 'MHO' => 'X',
			 'MLZ' => 'X',
			 'NAD' => 'X',
			 'NAL' => 'X',
			 'NH2' => 'X',
			 'NIT' => 'X',
			 'NLE' => 'X',
			 'ODS' => 'X',
			 'OXY' => 'X',
			 'PHD' => 'X',
			 'PHL' => 'X',
			 'PNL' => 'X',
			 'PPH' => 'X',
			 'PPL' => 'X',
			 'PRN' => 'X',
			 'PSS' => 'X',
			 'PSU' => 'X',
			 'PVL' => 'X',
			 'PY2' => 'X',
			 'QND' => 'X',
			 'QUO' => 'X',
			 'SEC' => 'X',
			 'SEG' => 'X',
			 'SEM' => 'X',
			 'SET' => 'X',
			 'SIN' => 'X',
			 'SLE' => 'X',
			 'THC' => 'X',
			 'TPN' => 'X',
			 'TRF' => 'X',
			 'UNK' => 'X',
			 'VAS' => 'X',
			 'YRR' => 'X',
			);

    my %fullname_to_one = ( 'GLYCINE'          => 'G',
			    'ALANINE'          => 'A',
			    'VALINE'           => 'V',
			    'LEUCINE'          => 'L',
			    'ISOLEUCINE'       => 'I',
			    'PROLINE'          => 'P',
			    'CYSTEINE'         => 'C',
			    'METHIONINE'       => 'M',
			    'HISTIDINE'        => 'H',
			    'PHENYLALANINE'    => 'F',
			    'TYROSINE'         => 'Y',
			    'TRYPTOPHAN'       => 'W',
			    'ASPARAGINE'       => 'N',
			    'GLUTAMINE'        => 'Q',
			    'SERINE'           => 'S',
			    'THREONINE'        => 'T',
			    'LYSINE'           => 'K',
			    'ARGININE'         => 'R',
			    'ASPARTATE'        => 'D',
			    'GLUTAMATE'        => 'E',
			    'ASPARTIC ACID'    => 'D',
			    'GLUTAMATIC ACID'  => 'E',
			    'ASPARTIC_ACID'    => 'D',
			    'GLUTAMATIC_ACID'  => 'E',
			    'SELENOMETHIONINE' => 'M',
			    'SELENOCYSTEINE'   => 'M',
			    'ADENINE'          => '0',
			    'CYTOSINE'         => '1',
			    'GUANINE'          => '2',
			    'THYMINE'          => '3',
			    'URACIL'           => '4'
			  );
    
    # map it
    #
    if (length $incode == 1) {
	$newcode = $one_to_three{$incode};
    }
    elsif (length $incode == 3) {
        $newcode = $three_to_one{$incode};
    }
    else {
	$newcode = $fullname_to_one{$incode};
    }


    # check for weirdness
    #
    if (! defined $newcode) {
#	&abort ("unknown residue '$incode'");
#	print STDERR ("unknown residue '$incode' (mapping to 'Z')\n");
#	$newcode = 'Z';
	if (! $silent) {
	    print STDERR ("unknown residue '$incode' (mapping to 'X')\n");
	}
	$newcode = 'X';
    }
    elsif ($newcode eq 'X') {
	if (! $silent) {
	    print STDERR ("strange residue '$incode' (seen code, mapping to 'X')\n");
	}
    }

    return $newcode;
}		      


# resAtomRecs()
#
sub resAtomRecs {
    my $code1 = shift;

    my %atomRecs = ( 
'A' => q{
ATOM   atmi  N   ALA  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  ALA  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   ALA  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   ALA  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  ALA  resi       0.000   0.000   0.000 -1.00  0.00
},
'C' => q{
ATOM   atmi  N   CYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  CYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   CYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   CYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  CYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  SG  CYS  resi       0.000   0.000   0.000 -1.00  0.00
},
'D' => q{
ATOM   atmi  N   ASP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  ASP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   ASP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   ASP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  ASP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  ASP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OD1 ASP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OD2 ASP  resi       0.000   0.000   0.000 -1.00  0.00
},
'E' => q{
ATOM   atmi  N   GLU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  GLU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   GLU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   GLU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  GLU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  GLU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD  GLU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OE1 GLU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OE2 GLU  resi       0.000   0.000   0.000 -1.00  0.00
},
'F' => q{
ATOM   atmi  N   PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD1 PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD2 PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE1 PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE2 PHE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CZ  PHE  resi       0.000   0.000   0.000 -1.00  0.00
},
'G' => q{
ATOM   atmi  N   GLY  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  GLY  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   GLY  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   GLY  resi       0.000   0.000   0.000 -1.00  0.00
},
'H' => q{
ATOM   atmi  N   HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  ND1 HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD2 HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE1 HIS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  NE2 HIS  resi       0.000   0.000   0.000 -1.00  0.00
},
'I' => q{
ATOM   atmi  N   ILE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  ILE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   ILE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   ILE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  ILE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG1 ILE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG2 ILE  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD1 ILE  resi       0.000   0.000   0.000 -1.00  0.00
},
'K' => q{
ATOM   atmi  N   LYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  LYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   LYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   LYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  LYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  LYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD  LYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE  LYS  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  NZ  LYS  resi       0.000   0.000   0.000 -1.00  0.00
},
'L' => q{
ATOM   atmi  N   LEU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  LEU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   LEU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   LEU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  LEU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  LEU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD1 LEU  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD2 LEU  resi       0.000   0.000   0.000 -1.00  0.00
},
'M' => q{
ATOM   atmi  N   MET  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  MET  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   MET  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   MET  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  MET  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  MET  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  SD  MET  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE  MET  resi       0.000   0.000   0.000 -1.00  0.00
},
'N' => q{
ATOM   atmi  N   ASN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  ASN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   ASN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   ASN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  ASN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  ASN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OD1 ASN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  ND2 ASN  resi       0.000   0.000   0.000 -1.00  0.00
},
'P' => q{
ATOM   atmi  N   PRO  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  PRO  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   PRO  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   PRO  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  PRO  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  PRO  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD  PRO  resi       0.000   0.000   0.000 -1.00  0.00
},
'Q' => q{
ATOM   atmi  N   GLN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  GLN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   GLN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   GLN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  GLN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  GLN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD  GLN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OE1 GLN  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  NE2 GLN  resi       0.000   0.000   0.000 -1.00  0.00
},
'R' => q{
ATOM   atmi  N   ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD  ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  NE  ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CZ  ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  NH1 ARG  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  NH2 ARG  resi       0.000   0.000   0.000 -1.00  0.00
},
'S' => q{
ATOM   atmi  N   SER  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  SER  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   SER  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   SER  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  SER  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OG  SER  resi       0.000   0.000   0.000 -1.00  0.00
},
'T' => q{
ATOM   atmi  N   THR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  THR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   THR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   THR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  THR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OG1 THR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG2 THR  resi       0.000   0.000   0.000 -1.00  0.00
},
'V' => q{
ATOM   atmi  N   VAL  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  VAL  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   VAL  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   VAL  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  VAL  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG1 VAL  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG2 VAL  resi       0.000   0.000   0.000 -1.00  0.00
},
'W' => q{
ATOM   atmi  N   TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD1 TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD2 TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  NE1 TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE2 TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE3 TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CZ2 TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CZ3 TRP  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CH2 TRP  resi       0.000   0.000   0.000 -1.00  0.00
},
'Y' => q{
ATOM   atmi  N   TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CG  TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD1 TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CD2 TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE1 TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CE2 TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CZ  TYR  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  OH  TYR  resi       0.000   0.000   0.000 -1.00  0.00
}
    );
    
    my $lines = $atomRecs{$code1};
    $lines =~ s/^\s+|\s+$//g;

    return split (/\n/, $lines);
}


# resAtomRecsNoSidechain()
#
sub resAtomRecsNoSidechain {
    my $code1 = shift;
    my %code3 = (
		 'A' => 'ALA',
		 'C' => 'CYS',
		 'D' => 'ASP',
		 'E' => 'GLU',
		 'F' => 'PHE',
		 'G' => 'GLY',
		 'H' => 'HIS',
		 'I' => 'ILE',
		 'K' => 'LYS',
		 'L' => 'LEU',
		 'M' => 'MET',
		 'N' => 'ASN',
		 'P' => 'PRO',
		 'Q' => 'GLN',
		 'R' => 'ARG',
		 'S' => 'SER',
		 'T' => 'THR',
		 'V' => 'VAL',
		 'W' => 'TRP',
		 'Y' => 'TYR'
		);

    my $lines = qq{
ATOM   atmi  N   cod  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CA  cod  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  C   cod  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  O   cod  resi       0.000   0.000   0.000 -1.00  0.00
ATOM   atmi  CB  cod  resi       0.000   0.000   0.000 -1.00  0.00
};
    $lines =~ s/cod/$code3{$code1}/g;
    $lines =~ s/^\s+|\s+$//g;

    return split (/\n/, $lines);
}


# getCommandLineOptions()
#
#  desc: get the command line options
#
#  args: none
#
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;
    local $usage = qq{
usage: $0 
\t -zonesfile       <zonesfile>
\t -fastafile       <fastafile>
\t -parentpdb       <parentpdb>
\t -outpdb          <outpdb>\
\t[-parentssa       <parentssa>]
\t[-outssa          <outssa>]
\t[-outcst          <outcst>]
\t[-takeoffpad      <T/F>]       (def: T)
\t[-loopregiononly  <T/F>]       (def: F)
\t[-sidechains      <T/F>]       (def: F)
\t[-keephetero      <T/F>]       (def: F)
};
#\t[-parentfssp   <parentfssp>]
#\t[-alignflag   <T/F>]     (def: T)

    # Get args
    #
    local %opts = ();
    &GetOptions (\%opts, "fastafile=s", 
		         "parentpdb=s", 
		         "parentssa=s", 
		         "zonesfile=s", 
		         "outpdb=s",
		         "outssa=s",
		         "outcst=s",
#		         "parentfssp=s",
		         "takeoffpad=s",
		         "loopregiononly=s",
		         "sidechains=s",
		         "keephetero=s");
		         #"alignflag=s", 

    # Check for legal invocation
    #
    if (! defined $opts{fastafile} ||
	! defined $opts{parentpdb} ||
	! defined $opts{zonesfile} ||
	! defined $opts{outpdb}
       ) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExistence ('f', $opts{fastafile});	
    &checkExistence ('f', $opts{parentpdb});	
    &checkExistence ('f', $opts{zonesfile});	#
    &checkExistence ('f', $opts{parentssa})  if ($opts{parentssa});	
#    &checkExistence ('f', $opts{parentfssp})  if ($opts{parentfssp});	

    # defaults
    #
    $opts{takeoffpad}     = 'T'    if (! defined $opts{takeoffpad});
    $opts{loopregiononly} = 'F'    if (! defined $opts{loopregiononly});
    $opts{sidechains}     = 'F'    if (! defined $opts{sidechains});
    $opts{keephetero}     = undef  if ($opts{keephetero} =~ /^F/i);
    #$opts{alignflag}      = 'T'    if (! defined $opts{alignflag});

    return %opts;
}

###############################################################################
# util
###############################################################################

sub maxInt {
    local ($v1, $v2) = @_;
    return ($v1 > $v2) ? $v1 : $v2;
}

sub tidyDecimals {
    my ($num, $decimal_places) = @_;
    if ($num !~ /\./) {
	$num .= '.' . '0' x $decimal_places;
	$num =~ s/^0+//;
    }
    else {
	if ($num =~ s/(.*\.\d{$decimal_places})(\d).*$/$1/) {
	    my $nextbit = $2;
	    if ($nextbit >= 5) {
		my $flip = '0.' . '0' x ($decimal_places - 1) . '1'; 
		$num += $flip;
	    }
        }
	$num =~ s/^0//;
	my $extra_places = ($decimal_places + 1) - length $num;
	$num .= '0' x $extra_places  if ($extra_places > 0);
    }

    return $num;
}

sub distsq {
    local @dims = @_;
    local $v = 0;
    foreach $dim (@dims) {
	$v += $dim*$dim;
    }
    return $v;
}

sub logMsg {
    local ($msg, $logfile) = @_;

    if ($logfile) {
        open (LOGFILE, ">".$logfile);
        select (LOGFILE);
    }
    else {
	select (STDERR);
    }
    print $msg, "\n";
    if ($logfile) {
        close (LOGFILE);
    }
    select (STDOUT);

    return 'true';
}

sub checkExistence {
    local ($type, $path) = @_;
    if ($type eq 'd') {
	if (! -d $path) { 
            print STDERR "$0: dirnotfound: $path\n";
            exit -3;
	}
    }
    elsif ($type eq 'f') {
	if (! -f $path) {
            print STDERR "$0: filenotfound: $path\n";
            exit -3;
	}
    }
}

sub abort {
    local $msg = shift;
    print STDERR "$0: $msg\n";
    exit -2;
}

sub writeBufToFile {
    ($file, $bufptr) = @_;
    if (! open (FILE, '>'.$file)) {
	&abort ("$0: unable to open file $file for writing");
    }
    print FILE join ("\n", @{$bufptr}), "\n";
    close (FILE);
    return;
}

sub fileBufString {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("$0: unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    local $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    return $buf;
}

sub fileBufArray {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("$0: unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    local $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    @buf = split (/$oldsep/, $buf);
    pop (@buf)  if ($buf[$#buf] eq '');
    return @buf;
}

###############################################################################
1; # package end
# end
###############################################################################
