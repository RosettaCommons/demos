#!/usr/bin/perl -w
##
##
## Copyright 2000, University of Washington
##   This document contains private and confidential information and
##   its disclosure does not constitute publication.  All rights are
##   reserved by University of Washington, except those specifically
##   granted by license.
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

#$width = 50;
#$spacing = 10;

###############################################################################
# init
###############################################################################

$| = 1;                                                   # don't buffer stdout

local %opts = &getCommandLineOptions ();
local $pdbfile = $opts{pdbfile};
local $chainOI = $opts{chain};

$pdbID = $pdbfile;
$pdbID =~ s!^.*?/?p?d?b?([^\.\/]+)\.[pe][dn][bt]\.?g?z?Z?$!$1!;
$pdbID =~ s/^(....)(.?)(.?).*/(lc$1).(uc$2).(uc$3)/e;

#$base_id = lc substr ($pdbID, 0, 4);

$chainOI = uc $chainOI;

###############################################################################
# main
###############################################################################

# read
#
@buf = &fileBufArray ($pdbfile);

$last_res_num  = -10000;
$chainOI_found = undef;
$start_res_num = undef;
$stop_res_num  = undef;
for ($i=0; $i <= $#buf; ++$i) {
    last if ($chainOI_found &&
	     ($buf[$i] =~ /^TER/ || $buf[$i] =~ /^ENDMDL/ || $buf[$i] =~ /^MODEL/));
    if ($buf[$i] =~ /^ATOM/ || $buf[$i] =~ /^HETATM/) {
	$chain   = substr ($buf[$i], 21, 1);
	$chain   = ($chain eq ' ') ? '_' : uc $chain;
	if (($chain eq '_' && $chainOI eq 'A') ||
	    ($chain eq 'A' && $chainOI eq '_')) {
	    print STDERR "$0: WARNING: changing chain sought from $chainOI to $chain\n";
	    $chainOI = $chain;
	}
	next if ($chain ne $chainOI);
	$chainOI_found = 'TRUE';
	$atom_type = substr ($buf[$i], 12, 4);
	if ($atom_type eq ' CA ') {
	    $res_num = substr ($buf[$i], 22, 5);     # includes insertion code
	    $res_num =~ s/\s+//g;
	    if ($res_num ne $last_res_num) {
		$last_res_num = $res_num;
		if (! defined $start_res_num) {
		    $start_res_num = $res_num;
		}
		$stop_res_num = $res_num;
		push (@fasta, &mapResCode (substr($buf[$i], 17, 3)));
	    }
	}
    }
}


# build fasta output
#
$outbuf = '';
#$outbuf .= ">$base_id$chainOI $start_res_num-$stop_res_num\n";
for ($i=0; $i <= $#fasta; ++$i) {
#    if ($i % $width == 0 && $i != 0) {
#	$outbuf .= "\n";
#    }
#    elsif ($i % $spacing == 0 && $i != 0) {
#	$outbuf .= ' ';
#    }
    $outbuf .= $i+1;
    $outbuf .= " ";
    $outbuf .= $fasta[$i];
    $outbuf .= " ";
    $outbuf .= "\.\n";
}
$outbuf =~ s/\s+$//g;
$outbuf .= "\n";


# output
print $outbuf;

# exit
#
exit 0;

###############################################################################
# subs
###############################################################################

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
    local $usage = qq{usage: $0 -pdbfile <pdbfile> [-chain <chain>]};

    # Get args
    #
    local %opts = ();
    &GetOptions (\%opts, "pdbfile=s", "chain=s");


    # Check for legal invocation
    #
    if (! defined $opts{pdbfile}) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExistence ('f', $opts{pdbfile});

    # defaults
    $opts{chain} = '_'  if (! defined $opts{chain});

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
    print $msg, "\n";
    if ($logfile) {
        close (LOGFILE);
        select (STDOUT);
    }
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
