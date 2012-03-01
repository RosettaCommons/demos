package hhsearch;

use alignment;

use constant MIN_E_VALUE => 1e-45; # must be greater than 0

sub parse_alignment_header {
  my $file = shift;

  open FILE, "<$file" or die $!;
  my @file = <FILE>;
  close FILE or die $!;

  my @columns = qw/ Prob E-value P-value Score SS Cols QueryHMM TemplateHMM TemplateLength/;
  my $in_alignment = 0;
  my @alignments;

  foreach my $line (@file) {
    chomp $line;
    if ( $line =~ /^\s*$/ ) {
      $in_alignment = 0;
      next; # skip blank lines
    }

    if ( $in_alignment ) {
      my $rest = substr($line,36);

      my @elements = split /\s+/, $rest;
      my %a;
      foreach my $i ( 0 .. scalar(@columns) - 1 ) {
        $a{ $columns[$i] } = $elements[$i];
      }
      push @alignments, \%a;
    }

    if ( $line =~ /No Hit/ ) {
      $in_alignment = 1;
    }
  }

  return \@alignments;
}

sub parse_alignments {
	my $file = shift;

	open FILE, "<$file" or die $!;
	my @file = <FILE>;
	close FILE or die $!;

	#print "reading alignments from $file!\n";

	my @alignments;
	my $aln;
	#print "read ", scalar(@file), " lines from $file.\n";
	foreach my $line (@file) {
		chomp $line;
		if ( $line =~ /^\s*$/ ) { next; }

		if ( $line =~ /^No (\d+)/ ) {
			if ( exists $aln->{rank} ) { push @alignments, $aln; }

			# make a new $aln
			my %nested_aln;
			$aln = \%nested_aln;
			$aln->{rank} = $1;
		} elsif ( $line =~ /Probab=([\d\.]+)\s+E-value=([+-e\d\.]+)\s+Score=([-\d\.]+)/ ) {
			$aln->{prob}    = $1;
			$aln->{e_value} = $2;
			$aln->{score}	 = $3;

			if ( $aln->{e_value} == 0 ) {
				$aln->{e_value} = MIN_E_VALUE;
			}
#Probab=100.00	E-value=0	Score=384.97	Aligned_columns=112	Identities=100%
		} elsif ( $line =~ /([QT]) (ss_pred)\s+([\w\-]+)/ ) {
			my $qt   = $1;
			my $type = $2;
			my $data = $3;
			$aln->{$qt}->{$type} .= $data;
		} elsif ( $line =~ /([QT]) (ss_dssp)\s+([\w\-]+)/ ) {
			my $qt   = $1;
			my $type = $2;
			my $data = $3;
			$aln->{$qt}->{$type} .= $data;
		} elsif ( $line =~ /([QT]) (ss_conf)\s+([\w\-]+)/ ) {
			my $qt   = $1;
			my $type = $2;
			my $data = $3;
			$aln->{$qt}->{$type} .= $data;
		} elsif ( $line =~ /([QT]) (Consensus)\s+([\w\-\~]+)/ ) {
			# do nothing for now
		#} elsif ( $line =~ /^([QT]) ([\(\)\d\w\.|_:\s]+)\s+(\d+)\s+([\w\-]+)\s+(\d+) \((\d+)\)/ ) {
		} elsif ( $line =~ /^([QT]) (.+)\s+(\d+)\s+([\w\-]+)\s+(\d+) \((\d+)\)/ ) {
			#Q ref|NP_244848.		1 MSFIEKMIGSLNDKREWKAMEARAKALPKEYHHAYKAIQKYMWTSGGPTDWQDTKRIFGGILDLFEEGAAEGKKVTDLTG	 80 (112)
			my $qt      = $1;
			my $name	   = $2;
			my $start   = $3;
			my $data    = $4;
			my $stop    = $5;
			my $length  = $6;

			#print "alignment = $data";

			$aln->{$qt}->{alignment} .= $data;
			# define start and stop appropriately
			if ( ! exists $aln->{$qt}->{start} ) {
				$aln->{$qt}->{start} = $start;
			}
			$aln->{$qt}->{stop} = $stop;
			$aln->{$qt}->{len}= $length;

			# remove trailing whitespace
			$name =~ s/\s+$//g;
			if ( $qt eq 'T' ) {
				if ( length($name) == 4 ) { $name .= '__'; }
				if (length($name) == 6) {
				    my $pdb = substr($name,0,4);
				    my $chain = substr($name,5);
				    $aln->{template_name} = "$pdb$chain";
				}
				else{
				    $aln->{template_name} = $name;
				}
			} else {
				$aln->{query_name} = $file;
			}
		} else {
			#warn "Error: don't recognize line:\n $line\n";
		}
	}

	push @alignments, $aln;

	my @objects = map { alignment->new( $_ ) } @alignments;
	return \@objects;
}

1;
