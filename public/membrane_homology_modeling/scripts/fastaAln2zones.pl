#!/usr/bin/perl
##
##  Initial Author: Bin QIan
##  $Revision: 1.2 $
##  $Date: 2004/01/02 22:34:39 $
##  $Author: bqian $
##
###############################################################################

###############################################################################
# conf
###############################################################################

$| = 1;                                              # disable stdout buffering

$debug = 0;                                                       # 0=off, 1=on

###############################################################################
# init
###############################################################################

# argv

if ($#ARGV != 1) {
    print STDERR "usage: $0 <fastafile alignment query / parent> <zonefile> \n";
    exit -1;
}

$alnfasta		= shift @ARGV;
$zonefile	 	= shift @ARGV;


###############################################################################
# main
###############################################################################

use Data::Dumper;

my @fasta = readFile($alnfasta);
my $query = undef;
$qaln='';
$saln='';
foreach (@fasta){
    chomp;
    if ($_ =~ />/){
        if ($query){
            $query = undef;
        }else{
            $query = 'true';
        }
    }else{
        if($query){
            $qaln .= $_;
        }else{
            $saln .= $_;
        }
    }
}

# get alignment between msa_fasta and db_fasta
#

#####
#clean up the alignment
@qalign_array = split(//, $qaln);
@salign_array = split(//, $saln);

$c = 0;
$qaln='';
$saln='';
foreach $q(@qalign_array){
  if ($q=~/\w/ || $salign_array[$c] =~/\w/){
	$qaln.=$q;
	$saln.=$salign_array[$c];
  }
  $c++;
}



my $qseq = $qaln;
my $sseq = $saln;
$qseq =~ s/\W//g;
$sseq =~ s/\W//g;
my $qend = length ($qseq);
my $send = length ($sseq);

my $l = length($qaln);
if ($l!=length($saln)) {print "aln length difference!\n";exit(1);}
my $qs = 1;
my $ss = 1;
my $qcount = 0;
my $scount = 0;
my $output  = "";
$inaln = 0;
$inaln = 1 if (substr($qaln,0,1)=~/\w/ && substr($saln,0,1)=~/\w/);
for (my $i=1; $i<=$l; $i++){
    if ((substr($qaln,$i-1,1) !~ /\w/ || substr($saln,$i-1,1) !~ /\w/) && $inaln){
	  $qe = $qcount;
      $se = $scount;
#      if ($qe < length ($qseq) && $se < length($sseq)) {
#         $qe = $qcount+$qstart-1;
#         $se = $scount+$sstart-1;
#      }
      if ($qe - $qs > 0){
        $output .= sprintf ("zone %4d-%-4d:%4d-%-4d\n",$qs,$qe,$ss,$se);
      }
	  $inaln = 0;
    }elsif ((substr($qaln,$i-1,1) =~ /\w/ && substr($saln,$i-1,1) =~ /\w/) && !$inaln){
	  $qs = $qcount+1;
      $ss = $scount+1;
#      if ($qs > 1 && $ss > 1) {
#        $qs = $qcount+$qstart+1;
#        $ss = $scount+$sstart+1;
#      }
	  $inaln = 1;
    }
    if (substr($qaln,$i-1,1) =~ /\w/){
	  $qcount++;
    }
    if (substr($saln,$i-1,1) =~ /\w/){
      $scount++;
    }
}
if ($inaln){
    $qe = $qend;
    $se = $send;
    if ($qe - $qs > 0){
        $output .= sprintf ("zone %4d-%-4d:%4d-%-4d\n",$qs,$qe,$ss,$se);
    }
}

open (Q, '>'.$zonefile);
print Q $output;

# done
exit 0;

###############################################################################
# subs
###############################################################################

# readFile
#
sub readFile {

    my $filename = shift;
    my @filecontent = ();
    if (!open(FILE,$filename)){&abort("fail open file $filename")}

    @filecontent = <FILE>;
    close(FILE);

    return @filecontent;
}


# make_blank_msa()
#
sub make_blank_msa {
    $out = '';
    @q_fasta = @{$db_fasta->{$query_id}};
    $seq = join ('', @{$db_fasta->{$query_id}});
    $range = '1-'.length($seq).':1-'.length($seq);
    $len_long_range = length ($range)  if (length ($range) > $len_long_range);
    $max_id_len = length ($query_id);

    # header
    $out .= "CODE".' 'x($max_id_len-length("CODE")).' ';
    $out .= sprintf ("%7s %5s %5s %5s %s ", 'LEN-ALN', 'IDENT', 'SCORE', 'E-VAL', (' 'x($len_long_range-length('RANGES'))).'RANGES');

    # query residue numbering
    $skips = 0;
    for ($i=0, $res_i=0; $i <= $#q_fasta; ++$i) {
	if ($q_fasta[$i] =~ /[a-zA-Z]/) {
	    ++$res_i;
	    if ($res_i % 10 == 0) {
		$out .= "|$res_i";
		$skips = length ($res_i);
	    } elsif ($skips != 0) {
		--$skips;
	    } else {
		$out .= ' ';
	    }
	} elsif ($skips != 0) {
	    --$skips;
	} else {
	    $out .= ' ';
	}
    }
    $out .= "\n";

    # data
    $out .= $query_id.' 'x($max_id_len-length($query_id));
    $out .= ' ';

    $len_aln = length ($seq);
    $identity = 100;
    $out .= sprintf ("%5s    %3s  ", $len_aln, $identity);
    $range = '1-'.length($seq).':1-'.length($seq);
    $out .=  sprintf ("%4s  %5s %s ", '***', '*****', (' 'x($len_long_range-length($range)).$range));
    $out .= "$seq\n";

    return $out;
}


# get completed msa
#
sub getMsaFullTrimHomologs {
    my ($msa_mapping, $query_aligned, $msa_fasta, $db_fasta, $query_id, @id_order) = @_;
    my $msa_full        = +{};
    my $last_mapped_res = +{};
    my $last_mapped_pos = +{};
    my $i = 0;
    my $id;


    for ($q_i=0; $q_i <= $#{$db_fasta->{$query_id}}; ++$q_i) {
	$msa_full->{$query_id}->[$q_i] = ($query_aligned->[$q_i])
	                                 ? $db_fasta->{$query_id}->[$q_i]
					 : lc $db_fasta->{$query_id}->[$q_i];
	foreach $id (@id_order) {
	    next if ($id eq $query_id);
	    next if (&listMember ($id, @bad_ids));
	    $mem_res_i = $msa_mapping->[$q_i]->{$id};
	    if (defined $mem_res_i) {
		$msa_full->{$id}->[$q_i] = $db_fasta->{$id}->[$mem_res_i];
		if (defined $last_mapped_res->{$id} &&
		    defined $last_mapped_pos->{$id} &&
		    $mem_res_i != $last_mapped_res->{$id}+1) {

		    $msa_full->{$id}->[$q_i] = lc $msa_full->{$id}->[$q_i];
		    $msa_full->{$id}->[$last_mapped_pos->{$id}] = lc $msa_full->{$id}->[$last_mapped_pos->{$id}];
		}
		$last_mapped_res->{$id} = $mem_res_i;
		$last_mapped_pos->{$id} = $q_i;
	    }
	    else {
		$msa_full->{$id}->[$q_i] = '.';
	    }
	}
    }

    # done!
    return $msa_full;
}


# get completed msa
#
sub getMsaFull {
    my ($msa_mapping, $query_aligned, $msa_fasta, $db_fasta, $query_id, @id_order) = @_;
    my $msa_full  = +{};
    my $i = 0;
    my $id;


    # build msa_full
    #
    my $last_res_i       = +{};
    my $last_mapped_pos  = +{} ;
    my $last_pos         = +{};
    for $id (@id_order) {
	next if (&listMember ($id, @bad_ids));
	$last_res_i->{$id} = -1;
	$last_mapped_pos->{$id} = -1;
    }
    for ($query_res_i=0, $i=0;
	 $query_res_i <= $#{$db_fasta->{$query_id}};
	 ++$query_res_i, ++$i) {

	# adjust to make space for insertions
	$max_diff = 0;
	foreach $id (@id_order) {
	    next if ($id eq $query_id);
	    next if (&listMember ($id, @bad_ids));
	    if (defined $msa_mapping->[$query_res_i]->{$id}) {
		$mem_res_i = $msa_mapping->[$query_res_i]->{$id};
		$space = $i - $last_mapped_pos->{$id} - 1;
		$num_missing_residues = $mem_res_i - $last_res_i->{$id} - 1;
		$diff = $num_missing_residues - $space;
		if ($last_res_i->{$id} != -1) {          # we'll do Nterm later
		    $max_diff = &maxInt ($max_diff, $diff);
		}
	    }
	}
	for ($j=0; $j < $max_diff; ++$j) {
	    foreach $id (@id_order) {
		next if (&listMember ($id, @bad_ids));
		$msa_full->{$id}->[$i+$j] = '.';
	    }
	}
	$i += $max_diff;

	# make insertions
	foreach $id (@id_order) {
	    next if ($id eq $query_id);
	    next if (&listMember ($id, @bad_ids));
	    if (defined $msa_mapping->[$query_res_i]->{$id}) {
		$mem_res_i = $msa_mapping->[$query_res_i]->{$id};
		if ($last_res_i->{$id} != -1
		    && $last_res_i->{$id} != $mem_res_i-1) {

		    for ($j=$i-1, $res_i=$mem_res_i-1;
			 $res_i > $last_res_i->{$id};
			 --$j, --$res_i) {

			$msa_full->{$id}->[$j] = lc $db_fasta->{$id}->[$res_i];
		    }
		}

		$last_res_i->{$id} = $mem_res_i;
		$last_mapped_pos->{$id} = $i;

		$mem_res = $db_fasta->{$id}->[$mem_res_i];
		$msa_full->{$id}->[$i] = $mem_res;
	    } else {
		$msa_full->{$id}->[$i] = '.';
	    }
	}

	$query_res = $db_fasta->{$query_id}->[$query_res_i];
	$msa_full->{$query_id}->[$i] = ($query_aligned->[$query_res_i])
	                                   ? $query_res
	                                   : lc $query_res;
	$last_pos->{$query_id} = $i;
    }


    # add Cterm
    $max_pos = $last_pos->{$query_id};
    foreach $id (@id_order) {
	next if ($id eq $query_id);
	next if (&listMember ($id, @bad_ids));
	$last_pos->{$id} = $last_mapped_pos->{$id};
	for ($mem_res_i=$last_res_i->{$id}+1, $i=$last_mapped_pos->{$id}+1;
	     $mem_res_i <= $#{$db_fasta->{$id}};
	     ++$mem_res_i, ++$i) {
	    $msa_full->{$id}->[$i] = lc $db_fasta->{$id}->[$mem_res_i];
	    $last_pos->{$id} = $i;
	    $max_pos = &maxInt ($max_pos, $last_pos->{$id});
	}
    }


    # fill in blanks at end
    foreach $id (@id_order) {
	next if (&listMember ($id, @bad_ids));
	for ($i=$last_pos->{$id}+1; $i <= $max_pos; ++$i) {
	    $msa_full->{$id}->[$i] = '.';
	}
    }


    # add Nterm
    my $first_res_i      = +{};
    my $first_mapped_pos = +{} ;
    my $max_Nterm_diff   = 0;
    foreach $id (@id_order) {
	next if ($id eq $query_id);
	next if (&listMember ($id, @bad_ids));
	$query_res_i = -1;
	for ($i=0; $i <= $#{$msa_full->{$id}}; ++$i) {
	    ++$query_res_i  if ($msa_full->{$query_id}->[$i] =~ /[a-zA-Z]/);
	    if (defined $msa_mapping->[$query_res_i]->{$id}) {
		$mem_res_i = $msa_mapping->[$query_res_i]->{$id};
		$space = $i;
		$num_missing_residues = $mem_res_i;
		$diff = $num_missing_residues - $space;
		$first_res_i->{$id} = $mem_res_i;
		$first_mapped_pos->{$id} = $i;
		$max_Nterm_diff = &maxInt ($max_Nterm_diff, $diff);
		last;
	    }
	}
    }
    foreach $id (@id_order) {
	next if (&listMember ($id, @bad_ids));
	for ($j=0; $j < $max_Nterm_diff; ++$j) {
	    unshift (@{$msa_full->{$id}}, '.');
	}
    }
    foreach $id (@id_order) {
	next if ($id eq $query_id);
	next if (&listMember ($id, @bad_ids));
	for ($mem_res_i=$first_res_i->{$id} - 1,
	     $i=$first_mapped_pos->{$id} + $max_Nterm_diff - 1;
	     $mem_res_i >= 0;
	     --$mem_res_i, --$i) {
	    $msa_full->{$id}->[$i] = lc $db_fasta->{$id}->[$mem_res_i];
	}
    }


    # done!
    return $msa_full;
}


# trimHomologs ()
#
sub trimHomologs {
    my ($msa_full, $query_id, @id_order) = @_;
    my $msa_new = +{};
    my $i, $j, $new_i;

    for ($i=0, $new_i=0; $i <= $#{$msa_full->{$query_id}}; ++$i) {
	if ($msa_full->{$query_id}->[$i] =~ /[A-Z]/i) {
	    foreach $id (@id_order) {
		next if (&listMember ($id, @bad_ids));
		$msa_new->{$id}->[$new_i] = $msa_full->{$id}->[$i];
	    }
	    ++$new_i;
	}
	else {
	    # lower case left and right side of insertion
	    foreach $id (@id_order) {
		next if (&listMember ($id, @bad_ids));
		if ($msa_full->{$id}->[$i] =~ /[a-z]/) {
		    for ($j=$i; $j >= 0; --$j) {
			if ($msa_full->{$id}->[$j] =~ /[A-Z]/) {
			    $msa_new->{$id}->[$new_i-1] = lc $msa_full->{$id}->[$j];
			    last;
			} else {
			    $msa_full->{$id}->[$j] = '.';
			}
		    }
		    for ($j=$i+1; $j <= $#{$msa_full->{$query_id}}; ++$j) {
			if ($msa_full->{$id}->[$j] =~ /[A-Z]/) {
			    $msa_full->{$id}->[$j] = lc $msa_full->{$id}->[$j];
			    last;
			} else {
			    $msa_full->{$id}->[$j] = '.';
			}
		    }
		}
	    }
	}
    }

    return $msa_new;
}


# align2seq()
#
sub align2seq {
    my ($scope, $seq1, $seq2, $par_id, $mem_id) = @_;
    my $alignment = +{};
    my $i;

    # straightforward linear ungapped (less expensive)
    #
    $alignment = &noGapLinearAlignment ($seq1, $seq2);

    if (! $alignment) {
	# gapping in seq2 alignment
	#
	$alignment = &alignSeqsNoGaps1 ($scope, $seq1, $seq2);

	# alignment must be perfect (that means > 1000*(n-1) residues)
	#
	if ($alignment->{aligned_residues_cnt} != $#{$seq2}+1) {
	    # view alignment
	    print STDERR "\n$mem_id:\n";

	    my $lastAlignment = -1;
	    print STDERR "seq1: ";
	    for ($i=0; $i <= $#{$seq1}; ++$i) {
		if (! defined $alignment->{$id}->{'1to2'}->[$i]
		    || $alignment->{$id}->{'1to2'}->[$i] == -1) {
		    print STDERR $seq1->[$i];
		} else {
		    print STDERR '-' while ($alignment->{$id}->{'1to2'}->[$i] > ++$lastAlignment);
		    print STDERR $seq1->[$i];
		}
	    }
	    print STDERR "\n";
	    $lastAlignment = -1;
	    print STDERR "seq2: ";
	    for ($i=0; $i <= $#{$seq2}; ++$i) {
		if (! defined $alignment->{$id}->{'2to1'}->[$i]
		    || $alignment->{$id}->{'2to1'}->[$i] == -1) {
		    print STDERR $seq2->[$i];
		} else {
		    print STDERR '-' while ($alignment->{$id}->{'2to1'}->[$i] > ++$lastAlignment);
		    print STDERR $seq2->[$i];
		}
	    }
	    print STDERR "\n";

	    # view scores
	    print STDERR "lenseq: ".($#{$seq2}+1)."\n";
	    print STDERR "score : ".$alignment->{aligned_residues_cnt}."\n";

	    print STDERR "$0: imperfect alignment between seq1 and seq2\n";
	}
    }

    return $alignment;
}


# noGapLinearAlignment ()
#
sub noGapLinearAlignment {
    my ($seq1, $seq2) = @_;
    my $alignment = undef;

    my $register_found = undef;
    my $register_shift = -1;
    my $seq1_copy = +[];

    for (my $i=0; $i <= $#{$seq1}; ++$i) {
	$seq1_copy->[$i] = $seq1->[$i];
    }

    for (my $base_i=0; $base_i <= $#{$seq1} - $#{$seq2} && ! $register_found; ++$base_i) {

	if (&getIdentity ($seq1_copy, $seq2) == 1) {
	    $register_found = 'true';
	    for (my $i=0; $i <= $#{$seq2}; ++$i) {
		$alignment->{'2to1'}->[$i] = $base_i + $i;
		$alignment->{'1to2'}->[$base_i+$i] = $i;
	    }
	    for (my $i=0; $i < $base_i; ++$i) {
		$alignment->{'1to2'}->[$i] = -1;
	    }
	    for (my $i=$base_i+$#{$seq2}+1; $i <= $#{$seq1}; ++$i) {
		$alignment->{'1to2'}->[$i] = -1;
	    }
	}
	else {
	    shift @{$seq1_copy};
	}
    }

    return $alignment;
}


# AlignSeqsNoGaps1()
#
sub alignSeqsNoGaps1 {
    my ($scope, $seq1, $seq2) = @_;
    my $alignment = +{};
    my $V = +[];                       # score for opt Q1...Qi<>T1...Tj
    #NOGAP my $E = +[];                # score for opt Q1...Qi<>T1...Tj, - <>Tj
    my $E = +[];                       # score for opt Q1...Qi<>T1...Tj, - <>Tj
    my $F = +[];                       # score for opt Q1...Qi<>T1...Tj, Qi<> -
    my $G = +[];                       # score for opt Q1...Qi<>T1...Tj, Qi<>Tj
    my $Vsource = +[];
# debug
#    my $pair = 1000;
#    my $mispair = -1000;
#    my $gap_init = 11;
#    my $gap_ext = 1;
    my $pair = 1000;
    my $mispair = -10000;
    my $gap_init_q = 750;
    #my $gap_init = 500;
    my $gap_init = 100;
    my $gap_ext = 0;
    my $INT_MIN = -100000;
    my ($i, $j);

# debug
#print "seq1: '";
#for ($i=0; $i<=$#{$seq1}; ++$i) {
#    print $seq1->[$i];
#}
#print "'\n";
#print "seq2: '";
#for ($i=0; $i<=$#{$seq2}; ++$i) {
#    print $seq2->[$i];
#}
#print "'\n";
#exit 0;
# end debug

    # basis
    #
    $V->[0]->[0] = 0;
    #$E->[0]->[0] = $INT_MIN;   # never actually accessed
    #$F->[0]->[0] = $INT_MIN;   # never actually accessed
    for ($i=1; $i <= $#{$seq1}+1; ++$i) {
	$V->[$i]->[0] = ($scope eq 'G')  ? -$gap_init - $i*$gap_ext  : 0;
	#NOGAP $E->[$i]->[0] = $INT_MIN;
	$E->[$i]->[0] = $INT_MIN;
    }
    for ($j=1; $j <= $#{$seq2}+1; ++$j) {
	$V->[0]->[$j] = ($scope eq 'G')  ? -$gap_init - $j*$gap_ext  : 0;
	$F->[0]->[$j] = $INT_MIN;
    }


    # recurrence
    #
    $maxScore = 0;
    for ($i=1; $i<= $#{$seq1}+1; ++$i) {
	for ($j=1; $j<= $#{$seq2}+1; ++$j) {
	    # note: seq1[i-1]==Qi and seq2[j-1]==Tj (i.e. res 1 stored as 0)

	    # G, F, and E
	    #
	    $G->[$i]->[$j] = $V->[$i-1][$j-1]
		+ &scorePair ($seq1->[$i-1], $seq2->[$j-1], $pair, $mispair);
	    $F->[$i]->[$j] = &maxInt ($V->[$i-1]->[$j] -$gap_init -$gap_ext, $F->[$i-1]->[$j] -$gap_ext);
	    #NOGAP $E->[$i]->[$j] = &maxInt ($V->[$i][$j-1] -$gap_init -$gap_ext, $E->[$i]->[$j-1] -$gap_ext);
	    $E->[$i]->[$j] = &maxInt ($V->[$i][$j-1] -$gap_init_q -$gap_ext, $E->[$i]->[$j-1] -$gap_ext);


	    # V
	    #
	    # Local scope and null string superior
	    #NOGAP if ($scope eq 'L' && $F->[$i]->[$j] < 0 && $E->[$i]->[$j] < 0 && $G->[$i]->[$j] < 0) {
	    if ($scope eq 'L' && $F->[$i]->[$j] < 0 && $E->[$i]->[$j] < 0 && $G->[$i]->[$j] < 0) {
	    #if ($scope eq 'L' && $F->[$i]->[$j] < 0 && $G->[$i]->[$j] < 0) {
		$V->[$i]->[$j] = 0;
		$Vsource->[$i]->[$j] = 'N';
	    }
	    # Global scope or null string inferior
	    else {
		#NOGAP if ($F->[$i]->[$j] >= $G->[$i]->[$j] || $E->[$i]->[$j] >= $G->[$i]->[$j]) {
		if ($F->[$i]->[$j] >= $G->[$i]->[$j] || $E->[$i]->[$j] >= $G->[$i]->[$j]) {
#		if ($F->[$i]->[$j] >= $G->[$i]->[$j]) {
		    #NOGAP if ($F->[$i]->[$j] >= $E->[$i]->[$j]) {
		    if ($F->[$i]->[$j] >= $E->[$i]->[$j]) {
			$V->[$i]->[$j] = $F->[$i]->[$j];
			$Vsource->[$i]->[$j] = 'F';
		    #NOGAP } else {
			#NOGAP $V->[$i]->[$j] = $E->[$i]->[$j];
			#NOGAP $Vsource->[$i]->[$j] = 'E';
		    #NOGAP }
		    } else {
			$V->[$i]->[$j] = $E->[$i]->[$j];
			$Vsource->[$i]->[$j] = 'E';
		    }
		} else {
		    $V->[$i]->[$j] = $G->[$i]->[$j];
		    $Vsource->[$i]->[$j] = 'G';
		}
	    }

	    # maxScore
	    #
	    if ($V->[$i]->[$j] > $maxScore) {
		$maxScore = $V->[$i]->[$j];
		$maxScore_i = $i;
		$maxScore_j = $j;
	    }
	}
    }
    $alignment->{score} = ($scope eq 'G')
	                      ? $V->[$#{$seq1}+1]->[$#{$seq2}+1]
			      : $maxScore;

    # walk back
    #
    if ($scope eq 'G') {
	$i = $#{$seq1}+1;
	$j = $#{$seq2}+1;
    } else {
	$i = $maxScore_i;
	$j = $maxScore_j;
    }
    $walkback_done = undef;
    while (! $walkback_done) {
	# Global stop condition
	if ($scope eq 'G' && ($i == 0 || $j == 0)) {
	    $walkback_done = 'TRUE';
	    last;
	}
	# Local stop condition
	elsif ($V->[$i]->[$j] == 0) {   # if(Vsource[i][j]=='N'): we hit null
	    $walkback_done = 'TRUE';
	    last;
	}
	if ($Vsource->[$i]->[$j] eq 'G') {
	    $alignment->{'1to2'}->[$i-1] = $j-1;   # seq1[i-1]==Qi
	    $alignment->{'2to1'}->[$j-1] = $i-1;   # seq2[j-1]==Tj
	    ++$alignment->{aligned_residues_cnt};
	    --$i; --$j;
	}
	elsif ($Vsource->[$i]->[$j] eq 'F') {
	    $alignment->{'1to2'}->[$i-1] = -1;
	    --$i;
	}
	else {
	    $alignment->{'2to1'}->[$j-1] = -1;
	    --$j;
	}
    }


    return ($alignment);
}


# scorePair()
#
sub scorePair {
    local ($r1, $r2, $pair, $mispair) = @_;
    return $pair  if ($r1 eq $r2);
    return $mispair;
}


# getIdentity()
#
sub getIdentity {
    my ($a, $b) = @_;
    local $score = 0;
    local $align_len = 0;
    local $a_len = $#{$a} + 1;
    local $b_len = $#{$b} + 1;
    local $len = ($a_len <= $b_len) ? $a_len : $b_len;
    &abort ("attempt to compare a zero length segment!")  if ($len == 0);

    for (my $i=0; $i < $len; ++$i) {
        next      if (! $a->[$i] || ! $b->[$i]);
	next      if ($a->[$i] eq ' ' || $b->[$i] eq ' ');
        next      if ($a->[$i] eq '-' || $b->[$i] eq '-');
	next      if ($a->[$i] eq '.' || $b->[$i] eq '.');
        ++$score  if ($a->[$i] eq $b->[$i]);
        ++$align_len;
    }
    &abort ("attempt to compare unaligned segment!")  if ($align_len == 0);
    return $score / $align_len;
}

###############################################################################
# util
###############################################################################

# getCommandLineOptions()
#
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;
    my $usage = qq{usage: $0
\t -id               <queryid>
\t -fastafile        <fastafile>
\t -blastfile        <blastfile>
\t[-m                <m_value>]                 (def: 0) (blast format: 0 or 6)
\t[-completehomologs <completehomologs (T/F)>]  (def: F)
\t[-trimhomologs     <trimhomologs (T/F)>]      (def: T)
\t[-db               <sequence_db>]             (def: /scratch/shared/genomes/nr)
\t[-outfile          <outfile>]                 (def: STDOUT)
\t[-pdbrecs          <pdbrecs>]                 (note: onlyif m=0)
\t[-nexttolast       <detection_block (T/F)>]   (def: F)
\t[-reclimit         <reclimit>]                (def: none)
};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, "id=s", "fastafile=s", "blastfile=s", "m=i", "completehomologs=s", "trimhomologs=s", "db=s", "outfile=s", "pdbrecs=s", "nexttolast=s", "reclimit=i");

    # Check for legal invocation
    if (! defined $opts{id} ||
	! defined $opts{fastafile} ||
	! defined $opts{blastfile}
        ) {
        print STDERR "$usage\n";
        exit -1;
    }

    # defaults
    #
    $opts{m}  = '0'                           if (! $opts{m});
    $opts{db} = "/scratch/shared/genomes/nr"  if (! $opts{db});
    $opts{trimhomologs} = 'TRUE'              if (! $opts{trimhomologs});

    # remappings
    $opts{completehomologs} = undef  if ($opts{completehomologs} =~ /^F/i ||
					 $opts{completehomologs} =~ /^N/i);
    $opts{trimhomologs} = undef      if ($opts{trimhomologs} =~ /^F/i ||
					 $opts{trimhomologs} =~ /^N/i);
    $opts{nexttolast} = undef        if ($opts{nexttolast} =~ /^F/i ||
					 $opts{nexttolast} =~ /^N/i);

    # existence checks
    #
    &checkExist ('f', $opts{db})  if ($opts{completehomologs});
    &checkExist ('f', $opts{fastafile});
    &checkExist ('f', $opts{blastfile});

    # secondary checks
    &abort ("-m must be 0 or 6")  if ($opts{m} != 0 && $opts{m} != 6);
    &abort ("cannot both complete and trim homologs")  if ($opts{completehomologs} && $opts{trimhomologs});
    &abort ("-m must be 0 if want pdbrecs")  if ($opts{m} != 0 && $opts{pdbrecs});

    return %opts;
}


# maxInt()
#
sub maxInt {
    local ($v1, $v2) = @_;
    return ($v1 > $v2) ? $v1 : $v2;
}


# listMember()
#
sub listMember {
    my ($item, @list) = @_;
    my $element;
    foreach $element (@list) {
	return $item  if ($item eq $element);
    }
    return undef;
}


# logMsg()
#
sub logMsg {
    my ($msg, $logfile) = @_;

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

# checkExist()
#
sub checkExist {
    my ($type, $path) = @_;
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
	elsif (! -s $path) {
            print STDERR "$0: emptyfile: $path\n";
            exit -3;
	}
    }
}

# runCmd
#
sub runCmd {
    my $cmd = shift;
    my $date = `date +'%Y-%m-%d_%T'`;  chomp $date;
    print "[$date]:$0: $cmd\n" if ($debug);
    &abort ("failure running cmd: $cmd\n") if (system ($cmd) != 0);
    #if (system ($cmd) != 0) {
    #    print STDERR ("failure running cmd: $cmd\n");
    #    return -2;
    #}
    return 0;
}


# abort()
#
sub abort {
    my $msg = shift;
    print STDERR "$0: $msg\n";
    exit -2;
}

# writeBufToFile()
#
sub writeBufToFile {
    my ($file, $bufptr) = @_;
    if (! open (FILE, '>'.$file)) {
	&abort ("$0: unable to open file $file for writing");
    }
    print FILE join ("\n", @{$bufptr}), "\n";
    close (FILE);
    return;
}

# fileBufString()
#
sub fileBufString {
    my $file = shift;
    my $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("$0: unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    my $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    return $buf;
}

# fileBufArray()
#
sub fileBufArray {
    my $file = shift;
    my $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("$0: unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    my $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    @buf = split (/$oldsep/, $buf);
    pop (@buf)  if ($buf[$#buf] eq '');
    return @buf;
}

# bigFileBufArray()
#
sub bigFileBufArray {
    my $file = shift;
    my $buf = +[];
    if ($file =~ /\.gz|\.Z/) {
        if (! open (FILE, "gzip -dc $file |")) {
            &abort ("$0: unable to open file $file for gzip -dc");
        }
    }
    elsif (! open (FILE, $file)) {
        &abort ("$0: unable to open file $file for reading");
    }
    while (<FILE>) {
        chomp;
        push (@$buf, $_);
    }
    close (FILE);
    return $buf;
}

###############################################################################
# end
###############################################################################
