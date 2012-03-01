package alignment;

use blast;
use hhsearch;
use File::Basename;
use List::Util qw/ min max /;
use FindBin;

sub filt_string {
	my $self = shift;
	my $this_query_name    = $self->get_query_name;
	my $this_template_name = $self->get_template_name;
	my $query_start        = $self->get_query_start    - 1;
	my $template_start     = $self->get_template_start - 1;
	my $query_sequence     = $self->get_query_alignment;
	my $template_sequence  = $self->get_template_alignment;

	my $score = 0;
	if ( exists $self->{score} ) {
		$score = $self->{score};
	}
	if ( exists $self->{e_value} ) {
		$score .= ' ' . $self->get_e_value;
	}

	if ( !defined $this_query_name ) {
		$this_query_name = 't000_';
	}

	my $pct_id  = sprintf( "%2.2f", 100 * $self->get_percent_identity );
	my $comment = join ' ', ( 'filt alignment', "($pct_id% id)" );
	if ( exists $self->{comment} ) {
		$comment = $self->{comment};
	}

	my $string = <<ALIGN;
## $this_query_name $this_template_name
# $comment
scores_from_program: $score
$query_start $query_sequence
$template_start $template_sequence
--
ALIGN

	return $string;
}

sub ungapped_query_len {
	my $self = shift;
	my $seq  = $self->get_query_alignment;
	$seq =~ s/\-//g;
	return length($seq);
}

sub ungapped_template_len {
	my $self = shift;
	my $seq  = $self->get_template_alignment;
	$seq =~ s/\-//g;
	return length($seq);
}

sub tex_formatted_string {
	my $self = shift;
	my $this_query_name    = $self->get_query_name;
	my $this_template_name = $self->get_template_name;
	my $query_start        = $self->get_query_start;
	my $template_start     = $self->get_template_start;
	my $query_sequence     = $self->get_query_alignment;
	my $template_sequence  = $self->get_template_alignment;

	my $score = 0;
	if ( exists $self->{score} ) {
		$score = $self->{score};
	}
	if ( exists $self->{e_value} ) {
		$score .= ' ' . $self->get_e_value;
	}

	my $comment = 'tex-formatted alignment';

	my $string = <<ALIGN;
# $comment
$this_query_name $query_start $query_sequence
$this_template_name $template_start $template_sequence
--
ALIGN

	return $string;
}


sub parse_alignments {
  my $file = shift;

  open FILE, "<$file" or die "Error opening file $file ($!)";
  my @file = <FILE>;
  close FILE or die $!;

  my @alignments;
  my $aln;
  foreach my $line (@file) {
    chomp $line;
    if ( $line =~ /^\s*$/ ) { next; }

    if ( $line =~ /^##\s+(.*)\s+([\d\w_]+)/ ) {
      if ( exists $aln->{query_name} ) { push @alignments, $aln; }
      # make a new $aln
      my %nested_aln;
      $aln = \%nested_aln;
      $aln->{query_name}    = $1;
      $aln->{template_name} = $2;
    #} elsif ( $line =~ /^(.*):\s*([-\d\.]+)/ ) {
    #  my $type = $1;
    #  my $data = $2;
    #  $aln->{$type} = $data;
	 } elsif ( $line =~ /^#(.*)/ && $line !~ /^##/ ) {
		my $comment = $1;
		$comment =~ s/^\s+//g;
		$comment =~ s/\s+$//g;
		$aln->{comment} = $comment;
	 } elsif ( $line =~ /scores_from_program:\s+(.*)/ ) {
		my @d = split /\s+/, $1;

		if ( @d ) { $aln->{score}   = shift @d; }
		if ( @d ) { $aln->{e_value} = shift @d; }

    } elsif ( $line =~ /^\s*(\d+)\s+([\w-]+)$/ ) {
      #print $line, "\n";
      my $start = $1;
      my $seq   = $2;
      if ( !exists $aln->{Q} ) {
        $ungapped_seq = $seq;
        $ungapped_seq =~ s/-//g; # remove gaps and calculate lengths

        $aln->{Q}->{start}     = $start + 1;
        $aln->{Q}->{stop}      = $start + length($ungapped_seq);
        $aln->{Q}->{alignment} = $seq;
        #print "setting query alignment to $seq!\n";
      } else {
        $ungapped_seq = $seq;
        $ungapped_seq =~ s/-//g; # remove gaps and calculate lengths

        $aln->{T}->{start}     = $start + 1;
        $aln->{T}->{stop}      = $start + length($ungapped_seq);
        $aln->{T}->{alignment} = $seq;
        #print "setting template alignment to $seq!\n";
      }
    }
  } # foreach my $line (@file)

  push @alignments, $aln;

  my @objects = map { alignment->new( $_ ) } @alignments;
  return \@objects;
}

sub new {
  my $class = shift;
  my $self  = shift;

  if ( !defined $self ) {
    die "Error: no arguments provided for alignment::new!\n";
  }

  bless $self, $class;
}

sub parse_aln {
	my $fn   = shift;
	my $type = shift;

	my @alignments;
	if ( $type eq 'hhsearch' ) {
		@alignments = hhsearch::parse_alignments( $fn );
	} elsif ( $type eq 'blast' ) {
		@alignments = blast::parse_alignments( $fn );
	} else {
		die "Error: don't know about format $type!\n";
	}

	return \@alignments;
}

sub parse_alignments_super_general {
	my $fn = shift;

	open FILE, "<$fn" or die $!;
	my @file = <FILE>;
	close FILE or die $!;

	my @alignments;
	my $aln;
	foreach my $line (@file) {
		chomp $line;
		if ( $line =~ /^\s*$/ ) { next; } # skip blank lines

		if ( $line =~ /^--$/ ) {
			if ( exists $aln->{query_name} ) {
				my $q_align     = $aln->{Q}->{alignment};
				my $t_align     = $aln->{T}->{alignment};
				my $q_rep_count = 0;
				my $t_rep_count = 0;
				while ( $q_align =~ s/^-//g ) { $q_rep_count++; }
				while ( $t_align =~ s/^-//g ) { $t_rep_count++; }
				my $max_rep     = max( $q_rep_count, $t_rep_count );
				my $q_diff = $max_rep - $q_rep_count;
				my $t_diff = $max_rep - $t_rep_count;
				$aln->{Q}->{start} = $aln->{Q}->{start} + ( $max_rep - $q_rep_count );
				$aln->{T}->{start} = $aln->{T}->{start} + ( $max_rep - $t_rep_count );
				$aln->{Q}->{alignment} = substr( $q_align, $q_diff );
				$aln->{T}->{alignment} = substr( $t_align, $t_diff );

				push @alignments, $aln;
			}
			# make a new $aln
			my %nested_aln;
			$aln = \%nested_aln;
		} elsif ( $line =~ /\s*(.*)\s+(\d+)\s+([\w-]+)$/ ) {
			my $name  = $1;
			my $start = $2;
			my $seq   = $3;

			$name =~ s/^\s+//g;
			$name =~ s/\s+$//g;

			if ( !exists $aln->{Q} ) {
				$aln->{query_name} = $name;
				$ungapped_seq = $seq;
				$ungapped_seq =~ s/-//g; # remove gaps and calculate lengths

				$aln->{Q}->{start}     = $start;
				#$aln->{Q}->{stop}      = $start + length($ungapped_seq);
				$aln->{Q}->{alignment} = $seq;
				#print "setting query alignment to $seq!\n";
			} else {
				$aln->{template_name} = $name;
				$ungapped_seq = $seq;
				$ungapped_seq =~ s/-//g; # remove gaps and calculate lengths

				$aln->{T}->{start}     = $start;
				#$aln->{T}->{stop}      = $start + length($ungapped_seq);
				$aln->{T}->{alignment} = $seq;
				#print "setting template alignment to $seq!\n";
			}
		} elsif ( $line =~ /score:/ ) {
			$line =~ s/://g;
			$line =~ s/^\s+//g;
			my ($name,$value) = split /\s+/, $line;
			$aln->{$name} = $value;
		} else {
			warn "line doesn't belong: $line\n";
		}
	} # foreach my $line (@file)

	if ( $aln->{Q}->{alignment} ) {
		my $q_align     = $aln->{Q}->{alignment};
		my $t_align     = $aln->{T}->{alignment};
		my $q_rep_count = 0;
		my $t_rep_count = 0;
		while ( $q_align =~ s/^-//g ) { $q_rep_count++; }
		while ( $t_align =~ s/^-//g ) { $t_rep_count++; }
		my $max_rep     = max( $q_rep_count, $t_rep_count );
		my $q_diff = $max_rep - $q_rep_count;
		my $t_diff = $max_rep - $t_rep_count;
		$aln->{Q}->{start} = $aln->{Q}->{start} + ( $max_rep - $q_rep_count );
		$aln->{T}->{start} = $aln->{T}->{start} + ( $max_rep - $t_rep_count );
		$aln->{Q}->{alignment} = substr( $q_align, $q_diff );
		$aln->{T}->{alignment} = substr( $t_align, $t_diff );
		push @alignments, $aln;
	}

	my @objects = map { alignment->new( $_ ) } @alignments;
	map { $_->recalculate_align_stops } @alignments;
	return \@objects;
}

sub parse_alignment_general {
	my $fn = shift;

	open FILE, "<$fn" or die $!;
	my @file = <FILE>;
	close FILE or die $!;

	my @alignments;
	my $aln;
	foreach my $line (@file) {
		chomp $line;
		if ( $line =~ /^\s*$/ ) { next; }

		if ( $line =~ /^##\s+(.*)\s+([\d\w_]+)/ ) {
			if ( exists $aln->{query_name} ) { push @alignments, $aln; }
			# make a new $aln
			my %nested_aln;
			$aln = \%nested_aln;
			$aln->{query_name}    = $1;
			$aln->{template_name} = $2;
		} elsif ( $line =~ /^(.*):\s*([-\d\.]+)/ ) {
			my $type = $1;
			my $data = $2;
			$aln->{$type} = $data;
		} elsif ( $line =~ /^(\d+)\s+([\w-]+)$/ ) {
			#print $line, "\n";
			my $start = $1;
			my $seq   = $2;
			if ( !exists $aln->{Q} ) {
				$ungapped_seq = $seq;
				$ungapped_seq =~ s/-//g; # remove gaps and calculate lengths

				$aln->{Q}->{start}     = $start;
				$aln->{Q}->{stop}      = $start + length($ungapped_seq);
				$aln->{Q}->{alignment} = $seq;
				#print "setting query alignment to $seq!\n";
			} else {
				$ungapped_seq = $seq;
				$ungapped_seq =~ s/-//g; # remove gaps and calculate lengths

				$aln->{T}->{start}     = $start;
				$aln->{T}->{stop}      = $start + length($ungapped_seq);
				$aln->{T}->{alignment} = $seq;
				#print "setting template alignment to $seq!\n";
			}
		} elsif ( $line =~ /score:/ ) {
			my (undef,$name,$value) = split /\s+/, $line;
			$aln->{$name} = $value;
		}
	} # foreach my $line (@file)

	push @alignments, $aln;

	my @objects = map { alignment->new( $_ ) } @alignments;
	return \@objects;
}

# data access methods

sub e_value {
	my $self = shift;
	if ( @_ ) {
		$self->{e_value} = shift;
	}

	return $self->{e_value};
}

sub get_prob {
  my $self = shift;
  return $self->{prob};
}

sub get_e_value {
  my $self = shift;
  return $self->e_value;
}

sub get_evalue {
	my $self = shift;
	return $self->get_e_value;
}

sub length {
	my $self = shift;

	if ( !exists $self->{Q}->{alignment} ) {
		die "Error: no alignment defined!\n";
		#use Data::Dumper::Simple;
		#print Dumper( $self );
	}

	return length( $self->get_query_alignment );
}

sub ungapped_length {
	my $self = shift;

	if ( !exists $self->{Q}->{alignment} ) {
		die "Error: no alignment defined!\n";
		#use Data::Dumper::Simple;
		#print Dumper( $self );
	}

	return length( $self->get_ungapped_query_alignment );
}

sub set_e_value {
	my $self = shift;
	$self->e_value( $new_value );
}

sub score {
	my $self = shift;
	if ( @_ ) {
		$self->{score} = shift;
	}
	return $self->{score};
}

sub get_score {
  my $self = shift;
  return $self->{score};
}

sub get_rank {
  my $self = shift;
  return $self->{rank};
}

sub get_template_start {
	my $self = shift;
	return $self->template_start;
}

sub template_start {
  my $self = shift;
  if ( @_ ) {
		$self->{T}->{start} = shift;
  }
  return $self->{T}->{start};
}

sub get_template_stop {
  my $self = shift;
  return $self->{T}->{stop};
}

sub query_start {
	my $self = shift;
	if ( @_ ) {
		$self->{Q}->{start} = shift;
	}
	return $self->{Q}->{start};
}

sub get_query_start {
  my $self = shift;
  return $self->query_start;
}

sub get_query_stop {
  my $self = shift;
  return $self->{Q}->{stop};
}

sub get_template_alignment {
  my $self = shift;
  return $self->{T}->{alignment};
}

sub get_query_alignment {
  my $self = shift;
  return $self->{Q}->{alignment};
}

sub get_ungapped_query_alignment {
	my $self = shift;
	my $seq = $self->get_query_alignment;
	$seq =~ s/-//g;
	return $seq;
}

sub get_ungapped_template_alignment {
	my $self = shift;
	my $seq = $self->get_template_alignment;
	$seq =~ s/-//g;
	return $seq;
}

sub get_template_ss {
  my $self = shift;
  return $self->{Q}->{ss_pred};
}

sub get_query_ss {
  my $self = shift;
  return $self->{Q}->{ss_pred};
}

sub get_query_name {
  my $self = shift;
  return $self->{query_name};
}

sub get_template_name {
  my $self = shift;
  return $self->template_name;
}

sub query_name {
	my $self = shift;
	if ( @_ ) {
		$self->{query_name} = shift;
	}
  return $self->{query_name};
}

sub template_name {
	my $self = shift;
	if ( @_ ) {
		$self->{template_name} = shift;
	}
  return $self->{template_name};
}

sub get_query_length {
  my $self = shift;
  #return $self->{Q}->{len};
  return length $self->{Q}->{alignment};
}

sub get_template_length {
  my $self = shift;
  #return $self->{T}->{len};
  return length $self->{T}->{alignment};
}

sub get_query_total_length {
  my $self = shift;
  return $self->{Q}->{len};
}
sub get_template_total_length {
  my $self = shift;
  return $self->{T}->{len};
}

sub pad_alignments {
  my $self = shift;
  my $query = $self->get_query_alignment;
  my $templ = $self->get_template_alignment;

  my $leading_pad = max( $self->get_query_start, $self->get_template_start );
  my $trailing_pad = max(
    ($self->get_query_length - $self->get_query_stop ),
    ($self->get_template_length - $self->get_template_stop) );

  print "leading_pad = $leading_pad, trailing_pad = $trailing_pad\n";
  $query    = '-' x $leading_pad . $query . '-' x $trailing_pad;
  $template = '-' x $leading_pad . $templ . '-' x $trailing_pad;
  $self->{Q}->{alignment} = $query;
  $self->{T}->{alignment} = $template;
}

# get a map of template sequence indices to query sequence indices
sub map_template_query {
  my $self = shift;
  my %alignment_map;

  my $query = $self->get_query_alignment;
  my $templ = $self->get_template_alignment;

  my $q_index = $self->get_query_start;
  my $t_index = $self->get_template_start;

  # set all mappings from template to query to be gap by default
  for ( 1 .. $self->get_template_length ) {
    $alignment_map{$_} = '-';
  }

  my @q = split //, $query;
  my @t = split //, $templ;
  while ( @q && @t ) {
    my $qc = shift @q;
    my $tc = shift @t;
    if ( $tc ne '-' && $qc ne '-' ) {
      $alignment_map{$t_index} = $q_index;
    }

    if ( $qc ne '-' ) { $q_index++ };
    if ( $tc ne '-' ) { $t_index++ };
  }

  return \%alignment_map;
}

sub define_aa_map {
  my $self = shift;

  my $query = $self->get_query_alignment;
  my $templ = $self->get_template_alignment;

  my $q_index = $self->get_query_start;
  my $t_index = $self->get_template_start;

  my @q = split //, $query;
  my @t = split //, $templ;
  my %map;
  while ( @q && @t ) {
    my $qc = shift @q;
    my $tc = shift @t;
    if ( $tc ne '-' && $qc ne '-' ) {
    }

    $map{Q}{$q_index} = $qc;
    $map{T}{$t_index} = $tc;

    if ( $qc ne '-' ) { $q_index++ };
    if ( $tc ne '-' ) { $t_index++ };
  }
  $self->{aa_map} = \%map;
  return $self;
} # define_aa_map

sub get_template_aa {
  my $self = shift;
  my $resi = shift;

  if ( !defined $self->{aa_map} ) {
    $self->define_aa_map;
  }

  #print "resi = $resi, aa = ", $self->{aa_map}{T}{$resi}, "\n";

  return $self->{aa_map}{T}{$resi};
}


# map query to template
sub map_query_template {
  my $self = shift;
  my %alignment_map;

  my $query = $self->get_query_alignment;
  my $templ = $self->get_template_alignment;

  my $q_index = $self->get_query_start;
  my $t_index = $self->get_template_start;

  # set all mappings from template to query to be gap by default
  for ( 1 .. $self->get_template_length ) {
    $alignment_map{$_} = '-';
  }

  my @q = split //, $query;
  my @t = split //, $templ;
  while ( @q && @t ) {
    my $qc = shift @q;
    my $tc = shift @t;
    if ( $tc ne '-' && $qc ne '-' ) {
        $alignment_map{$q_index} = $t_index;
    }

    if ( $qc ne '-' ) { $q_index++ };
    if ( $tc ne '-' ) { $t_index++ };
  }

  return \%alignment_map;
}

# returns the distance to the nearest gap for a given residue index in an alignment. Assumes
# the gaps are denoted with - characters, and that $resi is a one-indexed position within
# the query alignment.
#
# example:
# my $dgap = get_dgap( $resi )
sub get_dgap {
  my $self = shift;
	my $resi = shift;

  if ( !defined $resi ) {
    die "Error: resi not defined!\n";
  }

  if ( !defined $self->{dgaps} ) { $self->define_dgaps; }

  if ( !defined $self->{dgaps}->[$resi] ) {
    warn "Error: not defined for $resi in get_dgap\n";
    exit 1;
  }
  return $self->{dgaps}->[$resi];
}

sub define_dgaps {
  my $self = shift;

  my @query_seq   = split //, $self->get_query_alignment;
  my @templ_seq   = split //, $self->get_template_alignment;
  my $query_start = $self->get_query_start;
  my $query_stop  = $self->get_query_stop;

  my @gap_positions;
  my $i = $query_start;
  while ( @query_seq && @templ_seq ) {
    my $q = shift @query_seq;
    my $t = shift @templ_seq;
    if ( $q eq '-' || $t eq '-' ) {
      push @gap_positions, $i;
    }

    if ( $q ne '-' ) { $i++ };
  }

  # deal with case of an ungapped alignment
  my @dgaps;
  if ( scalar(@gap_positions) == 0 ) {
    my $query_length = length( $self->get_query_alignment );
    $dgaps[$_] = $query_length for ( $query_start .. $query_stop );
  } else {
    for my $i ( $query_start .. $query_stop ) {
      my @sorted_positions = sort { abs( $i - $a ) <=> abs( $i - $b ) } @gap_positions;
      $dgaps[$i] = abs($i - $sorted_positions[0]);
    }
  }

  $self->{dgaps} = \@dgaps;
}

sub get_blosum_score {
	my $self = shift;
	my $resi = shift;

	if ( !defined $resi ) {
	  die "Error: resi not defined in get_blosum_score!\n";
	}

	if ( !defined $self->{blosum_scores} ) { $self->define_blosum_scores; }

	if ( !defined $self->{blosum_scores}->[$resi] ) {
	  return -4; # marker value for a gap.
	} else {
	  return $self->{blosum_scores}->[$resi];
	}
}

sub define_blosum_scores {
  my $self = shift;
  my $init = 0;
  open BLOSUM, "$FindBin::Bin/../cm_data/BLOSUM62" or die $!;
  my @b = <BLOSUM>;
  close BLOSUM or die $!;

  # skip comment lines
  @b = grep { !/^#/ } @b;
  my $header_line = shift @b;
  chomp $header_line;
  my @columns = grep { !/^\s*$/ } split /\s+/, $header_line;

  foreach my $line (@b) {
    chomp $line;
    my @elements = split /\s+/, $line;
    my $element_name = shift @elements;
    foreach my $i ( 0 .. scalar(@elements) - 1 ) {
      my $column_name = $columns[$i];
      my $element = $elements[$i];
      $blosum_data{$column_name}{$element_name} = $element;
    }
  }

  my @blosum_scores;
  my @query_seq   = split //, $self->get_query_alignment;
  my @templ_seq   = split //, $self->get_template_alignment;
  my $query_start = $self->get_query_start;
  my $query_stop  = $self->get_query_stop;

  my $i = $query_start;
  while ( @query_seq && @templ_seq ) {
    my $resi = shift @query_seq;
    my $resj = shift @templ_seq;
    if ( $resi ne '-' && $resj ne '-' ) {
      $blosum_scores[$i] = $blosum_data{$resi}{$resj};
      $i++;
    }
  }

  $self->{blosum_scores} = \@blosum_scores;
}

sub define_map_scores {
  my $self = shift;
  my $fn1  = shift;
  my $fn2  = shift;

  my $hhalign = './hhalign.new';

  my $atab_file = '/tmp/' . basename $fn1 . basename $fn2;
  my $cmd = "$hhalign -i $fn1 -t $fn2 -atab $atab_file -v 0 -glob";
  system( $cmd );

  open FILE, "<$atab_file" or die $!;
  my @lines = <FILE>;
  close FILE or die $!;

  my %map_scores;
  shift @lines;
  foreach my $line (@lines) {
    chomp $line;
    my (undef,$q_resi,$t_resi,$score,$ss,$probab) = split /\s+/, $line;
    $map_scores{ $q_resi } = $score;
    #$map_scores{ $q_resi } = $probab;
  }

  #system( "cp $atab_file ~tex/template_stats/hhsearch_runs/" );
  unlink $atab_file;
  return \%map_scores;
}

sub get_map_score {
  my $self = shift;
  my $resi = shift;
  return $self->{map_scores}[ $resi ];
}

sub get_gap_percentage {
  my $self = shift;
  my $query = $self->get_query_alignment;
  my $templ = $self->get_template_alignment;
  my $gaps = 0;
  for ($j = 0; $j < length($query); $j++) {
    if( substr($query,$j,1) eq "-" || substr($templ,$j,1) eq "-") {
      $gaps++;
    }
  }
  my $pergap = $gaps / length($query);
  return $pergap;
}

sub get_percent_identity {
  my $self = shift;

  my $query = $self->get_query_alignment;
  my $templ = $self->get_template_alignment;
  my $identities = 0;
  for ($j = 0; $j < length($query); $j++) {
    if( substr($query,$j,1) eq substr($templ,$j,1) ) {
      $identities++;
    }
  }

  my $percentage = $identities / length($query);
  return $percentage;
}

sub get_full_blosum_score {
  my $self = shift;

  $self->define_blosum_scores;
  my $score_tot = 0;
  my $count = 0;
  for my $score ( @{ $self->{blosum_scores} } ) {
    if ( defined $score ) {
      #print "score = $score!\n";
      $score_tot += $score;
      $count++;
    }
  }

  #return ($score_tot / $count);
  return $score_tot;
}

sub get_total_query_length {
	my $self = shift;
	if ( defined $self->{total_query_len} ) {
		return $self->{total_query_len};
	} else {
		die "Error: total_query_length not defined for alignment objects!";
	}
}

sub recalculate_align_stops {
	my $self = shift;
	my $ungapped_q = $self->get_query_alignment;
	my $ungapped_t = $self->get_template_alignment;
	$ungapped_q =~ s/-//g;
	$ungapped_t =~ s/-//g;

	$self->{Q}->{stop} = $self->get_query_start    + length($ungapped_q);
	$self->{T}->{stop} = $self->get_template_start + length($ungapped_t);
	return $self;
}

sub get_percent_query_aligned {
	my $self = shift;
	my $total_query_length = $self->get_total_query_length;
	my $query_seq = $self->get_query_alignment;
	$query_seq =~ s/-//g;

	return length($query_seq) / $total_query_length;
}

sub get_query_res {
	my $self = shift;
	my $resi = shift;

	my $idx = $resi - $self->get_query_start;
	if ( $idx <= $self->get_query_length ) {
		return substr( $self->get_query_alignment, $idx, 1 );
	}
	return '-';
}

sub get_template_res {
	my $self = shift;
	my $resi = shift;

	my $idx = $resi - $self->get_template_start;
	if ( $idx <= $self->get_template_length ) {
		return substr( $self->get_template_alignment, $idx, 1 );
	}
	return '-';
}

1;
