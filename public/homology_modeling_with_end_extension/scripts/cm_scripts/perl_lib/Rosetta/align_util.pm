package Rosetta::align_util;

use Exporter;
@ISA = qw/ Exporter /;

@EXPORT_OK = qw/
	convert_aln
	read_alns
	pairwise_percent_id
	align_two_seqs
	read_fasta
	read_fasta_as_hash
/;

sub convert_aln {
	use alignment;

	my $aln = shift;
	my $format_out = shift;

	my $str;
	if ( $format_out eq 'grishin' ) {
		$str = $aln->filt_string;
	} elsif ( $format_out eq 'tex' ) {
		$str = $aln->tex_formatted_string;
	} else {
		die "Error: don't recognize format $format!\n";
	}

	return $str;
}

sub pairwise_percent_id {
	my $seq1 = shift;
	my $seq2 = shift;

   use File::Temp qw/ tempfile /;
   use blast;

   my ($fh1, $fasta1)      = tempfile();
   my ($fh2, $fasta2)      = tempfile();
   my (undef, $blast_file) = tempfile();

	print $fh1 "> seq1\n$seq1\n";
	print $fh2 "> seq2\n$seq2\n";

   my $blast_cmd = "bl2seq -F f -i $fasta1 -j $fasta2 -p blastp -o $blast_file";
   system( $blast_cmd );

   if ( -z $blast_file ) {
      warn "Error: $blast_file has zero length!\n";
   }

   my $retval = 0;
   my @blast_aln = @{ blast::parse_alignments( $blast_file ) };
   if ( scalar(@blast_aln) != 0 ) {
      $retval = $blast_aln[0]->get_percent_identity;;
   } else {
		warn "Error: no alignments in $blast_file\n";
	}

   unlink( $fasta1 );
   unlink( $fasta2 );
   unlink( $blast_file );

   return $retval;
}

sub align_two_seqs {
	my $seq1 = shift;
	my $seq2 = shift;

   use File::Temp qw/ tempfile /;
   use blast;

   my ($fh1, $fasta1)      = tempfile();
   my ($fh2, $fasta2)      = tempfile();
   my (undef, $blast_file) = tempfile();

	print $fh1 "> seq1\n$seq1\n";
	print $fh2 "> seq2\n$seq2\n";

   my $blast_cmd = "bl2seq -F f -i $fasta1 -j $fasta2 -p blastp -o $blast_file";
   system( $blast_cmd );

   if ( -z $blast_file ) {
      warn "Error: $blast_file has zero length!\n";
   }

   my $retval = 0;
   my @blast_aln = @{ blast::parse_alignments( $blast_file ) };

   if ( scalar(@blast_aln) == 0 ) {
      #die "Error: can't align sequences!\n";
		return 0;
   }

   return $blast_aln[0];
}

sub parse_alignments_tex {
	my $fn = shift;

	my @alns;
	my $aln = {};
	open FILE, "<$fn" or die $!;
	while ( my $line = <FILE> ) {
		chomp $line;
		if ( $line =~ /query\s+(\d+)\s+([\w\-]+)/ ) {
			$aln->{Q}->{start} = $1;
			$aln->{Q}->{alignment} = $2;
		} elsif ( $line =~ /([\d\w\.]+)\s+(\d+)\s+([\w\-]+)/ ) {
			$aln->{template_name} = $1;
			$aln->{T}->{start} = $2;
			$aln->{T}->{alignment} = $3;
		} elsif ( $line =~ /^--$/ ) {
			if ( exists $aln->{Q}->{start} ) {
				push @alns, $aln;
			}
		} else {
			#print "don't recognize line:\n$line!\n";
		}
	}
	close FILE or die $!;
	if ( exists $aln->{Q}->{start} ) {
		push @alns, $aln;
	}

	my @objects = map { alignment->new($_) } @alns;

	return \@objects;
}

sub read_alns {
	my $aln_file    = shift;
	my $aln_format  = shift;
	my $ev_map_file = shift;

	my @alns;
	if ( $aln_format eq 'hhsearch' ) {
		@alns = @{ hhsearch::parse_alignments( $aln_file ) };
	} elsif ( $aln_format eq 'grishin' ) {
		@alns = @{ alignment::parse_alignments( $aln_file ) };
	} elsif ( $aln_format eq 'blast' ) {
		@alns = @{ blast::parse_alignments( $aln_file ) };
	} elsif ( $aln_format eq 'tex' ) {
		@alns = @{ parse_alignments_tex($aln_file) };
	} else {
		die "Error: don't recognize format $aln_format!\n";
	}

	if ( $ev_map_file && -f $ev_map_file ) {
		print STDERR "ev_map_file = $ev_map_file\n";
	}

	if ( $ev_map_file && -f $ev_map_file ) {
		open FILE, "<$ev_map_file" or die $!;
		my $header = <FILE>;
		my %ev_map;
		while ( my $line = <FILE> ) {
			chomp $line;
			my ($template,$e_value) = split /\s+/, $line;
			$ev_map{$template} = $e_value;
		}

		foreach my $aln (@alns) {
			my $template_name = join '', (
				lc( substr( $aln->get_template_name, 0, 4 ) ),
				uc( substr( $aln->get_template_name, 4, 1 ) )
			);
			$aln->template_name( $template_name );

			if ( exists $ev_map{$template_name} ) {
				#print "setting e-value of $template_name to ", $ev_map{$template_name}, "\n";
				$aln->e_value( $ev_map{$template_name} );
			}
		}
	}

	if ( wantarray ) {
		return @alns;
	}
	return \@alns;
}

sub read_fasta {
	my $fn = shift;

	my %seqs;
	my $current_id  = 0;
	my $current_seq = '';
	open FILE, "<$fn" or die $!;
	while ( my $line = <FILE> ) {
		chomp $line;
		if ( $line =~ /^>\s*(.*)/ ) {
			$current_id = $1;
		} else {
			$line =~ s/\s+//g;
			if ( $current_id ) {
				$seqs{$current_id} .= $line;
			}
		}
	}

	return \%seqs;
}

sub read_fasta_as_hash {
	my $fn = shift;

	my %seqs;
	my $current_id  = 0;
	my $current_seq = '';
	open FILE, "<$fn" or die $!;
	while ( my $line = <FILE> ) {
		chomp $line;
		if ( $line =~ /^>\s*(.*)/ ) {
			$current_id = $1;
		} else {
			$line =~ s/\s+//g;
			if ( $current_id ) {
				$seqs{$current_id} .= $line;
			}
		}
	}

	return \%seqs;
}

1;
