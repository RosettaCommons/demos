package CstUtil;

use Rosetta::PDB;
@ISA = qw/ Exporter /;

@EXPORT_OK = qw/
	dgaussian
	csts_to_str
	read_data
	data_to_hashref
	trim_whitespace
	read_csts
	merge_redundant_csts
	reweight_csts
	calc_probability
	calc_native_probability
	distances_from_pdb
	data_to_str
	data_from_str
	cst_set_to_str
/;

sub cst_set_to_str {
	my $csts = shift;

	my @lines;
	foreach my $resi ( keys %$csts ) {
	foreach my $resj ( keys %{ $csts->{$resi} } ) {
		my $line = "AtomPair CA $resi CA $resj ";
		my $n_csts = scalar( @{ $csts->{$resi}{$resj} } );
		$line .= "SOGFUNC $n_csts ";

		foreach my $cst ( @{ $csts->{$resi}{$resj} } ) {
			$line .= join ' ', 
				map { sprintf( "%.3f", $_ ) }
				( $cst->{mean}, $cst->{sdev}, $cst->{weight} );
			$line .= ' ';
		}
		push @lines, $line;
	}
	}

	my $str = join "\n", @lines;
	$str .= "\n";
	return $str;
}

sub calc_probability_prior {
	my $pdb_fn = shift;
	my $options = {};
	$options->{storage} = 0;

	if ( @_ ) {
		my $other_options = shift;
		foreach my $key ( keys %$other_options ) {
			$options->{$key} = $other_options->{$key};
		}
	}

	my $fasta = Rosetta::PDB::pdb_to_fasta( $pdb_fn );
	my @lines = grep { !/^>/ } split /\n/, $fasta;
	$fasta = join '', @lines;

	my %constraints;
	for my $resi ( 1 .. length($fasta) ) {
	for my $resj ( $resi + 1 .. length($fasta) ) {
		push @{ $constraints{$resi}{$resj} },
			{ mean => 18.3, sdev => 7.3, weight => 1.0 };
	}
	}

	return calc_probability( $pdb_fn, \%constraints, $options );
}

sub calc_native_probability {
	my $pdb_fn = shift;
	my $sdev   = shift;
	my $options = {};
	$options->{storage} = 0;

	if ( @_ ) {
		my $other_options = shift;
		foreach my $key ( keys %$other_options ) {
			$options->{$key} = $other_options->{$key};
		}
	}

	my $fasta = Rosetta::PDB::pdb_to_fasta( $pdb_fn );
	my @lines = grep { !/^>/ } split /\n/, $fasta;
	$fasta = join '', @lines;

	my %dists = %{ stored_distances( $options->{storage}, $pdb_fn ) };
	my %constraints;
	for my $resi ( 1 .. length($fasta) ) {
	for my $resj ( $resi + 1 .. length($fasta) ) {
		my $dist = $dists{$resi}{$resj};
		if ( $dist < 12 ) {
			push @{ $constraints{$resi}{$resj} }, 
				{ mean => $dist, sdev => $sdev, weight => 1.0 };
		} else {
			push @{ $constraints{$resi}{$resj} }, 
				{ mean => 18.3, sdev => 7.3, weight => 1.0 };
		}
	}
	}

	return calc_probability( $pdb_fn, \%constraints, $options );
}

sub maximize_p_native {
	my $dist = shift;
	my $csts = shift;

	push @$csts, { mean => 18.3, sdev => 7.3, weight => 0 };

	foreach my $cst (@$csts) {
		$cst->{prob} = dgaussian( $cst->{mean}, $cst->{sdev}, $dist );
	}

	my @new_csts = sort { $b->{prob} <=> $a->{prob} } @$csts;
	foreach my $idx ( 0 .. scalar(@new_csts) - 1 ) {
		if ( $idx == 0 ) {
			$new_csts[$idx]{weight} = 1;
		} else {
			$new_csts[$idx]{weight} = 0;
		}
		delete $new_csts[$idx]->{prob};
	}

	return \@new_csts;
}

sub stored_distances {
	my $storage = shift;
	my $pdb_fn  = shift;

	my %dists;
	if ( exists $storage->{$pdb_fn} ) {
		%dists = %{ $storage->{$pdb_fn} };
	} else {
		%dists = %{ distances_from_pdb( $pdb_fn ) };
		$storage->{$pdb_fn} = \%dists;
	}

	return \%dists;
}

sub calc_probability {
	my $pdb_fn    = shift;
	my $orig_csts = shift;

	use Storable qw/ dclone /;
	my $csts = dclone($orig_csts);

	my $options = {};
	$options->{storage} = 0;
	$options->{maximize_p_native} = 0;

	if ( @_ ) {
		my $other_options = shift;
		foreach my $key ( keys %$other_options ) {
			$options->{$key} = $other_options->{$key};
		}
	}

	my $fasta = Rosetta::PDB::pdb_to_fasta( $pdb_fn );
	my @lines = grep { !/^>/ } split /\n/, $fasta;
	$fasta = join '', @lines;

	my @pairs;
	for my $resi ( 1 .. length($fasta) ) {
	for my $resj ( $resi + 10 .. length($fasta) ) {
		push @pairs, [ $resi, $resj ];
	}
	}

	my %dists = %{ stored_distances( $options->{storage}, $pdb_fn ) };

	my $n_pairs  = 0;
	my $total_ll = 0;
	foreach my $pair (@pairs) {
		my $resi = $pair->[0];
		my $resj = $pair->[1];

		my $ll = 0;
		my $dist = $dists{$resi}{$resj};
		if ( !defined $dist ) {
			die "Error: distance between $resi,$resj not defined!\n";
		}

		if ( scalar( @{ $csts->{$resi}{$resj} } ) > 0 ) {
			my $weighted_prob = 0;
			my $total_weight  = 0;

			my @csts = @{ $csts->{$resi}{$resj} };

			if ( $options->{maximize_p_native} ) {
				@csts = @{ maximize_p_native( $dist, \@csts ) };
			}

			foreach my $cst (@csts) {
				if ( ($cst->{weight} - 1e-20) > 0 ) {
					$weighted_prob += $cst->{weight} *
						dgaussian( $cst->{mean}, $cst->{sdev}, $dist );
					$total_weight += $cst->{weight};
				}
			}

			if ( $weighted_prob > 1 || $weighted_prob <= 0 ) {
				use Data::Dumper::Simple;
				print Dumper( @csts );
				print "$resi, $resj\n";
				print "weights are : ", join ' ', map { $_->{weight} } @csts;
				die "total_weight = $total_weight\n";
			}
			#$ll = max( log( $weighted_prob ), - );
			$ll = log($weighted_prob);

			if ( abs( $total_weight - 1 ) > 1e-3 ) {
				print "$resi, $resj\n";
				print "weights are : ", join ' ', map { $_->{weight} } @csts;
				die "total_weight = $total_weight\n";
			}
		} else {
			my $bg_probability = log( dgaussian( 18.3, 7.3, $dist ) );
			$ll = $bg_probability;
		}

		$total_ll += $ll;
		$n_pairs++;
		#print "log_prob($resi,$resj,$dist) = $ll\n";
	}

	#print STDERR "n_pairs = $n_pairs\n";
	return $total_ll;
}

sub distances_from_pdb {
	my $pdb_fn = shift;

	use Rosetta;
	my $binary = "/work/tex/src/mini/bin/batch_distances";
	my $args = "-in:file:s $pdb_fn -mute all -database /work/tex/minirosetta_database";

	print STDERR "reading data from $pdb_fn ... ";
	my $mgr = Rosetta->new( compiler => 'gcc' );
	my $output = $mgr->run_program( $binary, $args );
	if ( -z $output ) {
		die "Error: no output from program!\n";
	}
	my @data = @{ data_from_str( $output ) };

	my %distances;
	foreach my $d (@data) {
		$distances{ $d->{resi_idx} }{ $d->{resj_idx} } = $d->{dist};
	}
	print STDERR "done.\n";
	return \%distances;
}

sub csts_to_str {
	my $resi = shift;
	my $resj = shift;
	my $csts = shift;

	my $str = "AtomPAIR CA $resi CA $resj SOGFUNC";
	foreach my $cst (@$csts) {
		$str = join ' ', ( $str, $cst->{mean}, $cst->{sdev}, $cst->{weight} );
	}

	return $str;
}

sub read_data {
	my $fn = shift;
	open FILE, "<$fn" or die "Error opening file $fn ($!)";
	my $header = <FILE>;
	$header = trim_whitespace( $header );
	my @columns = split /\s+/, $header;
	my @data;
	while ( my $line = <FILE> ) {
		$line = trim_whitespace( $line );
		my @d = split /\s+/, $line;

		push @data, data_to_hashref( \@d, \@columns );
	}
	close FILE or die $!;

	return \@data;
}

sub data_from_str {
	my $str = shift;

	my @lines = map { trim_whitespace($_) } split /\n/, $str;
	my $header = shift @lines;
	my @columns = split /\s+/, $header;

	my @data;
	foreach my $line (@lines) {
		my @d = split /\s+/, $line;
		push @data, data_to_hashref( \@d, \@columns );
	}
	return \@data;
}

sub data_to_str {
	my $data = shift;

	my @col_names = sort keys %{ $data->[0] };
	if ( @_ ) {
		my $ref = shift;
		@col_names = @$ref;
	}

	my $str = join ' ', @col_names;
	$str .= "\n";
	foreach my $d (@$data) {
		if ( grep { !exists $d->{$_} } @col_names ) {
			die "Error: mismatch between column names and data!\n";
		}
		$str .= join ' ', map { $d->{$_} } @col_names;
		$str .= "\n";
	}

	return $str;
}

sub data_to_hashref {
	my $d = shift;
	my $c = shift;

	if ( scalar(@$d) != scalar(@$c) ) {
		#use Data::Dumper::Simple;
		#print Dumper($d);
		#print Dumper($c);
		warn "columns are:\n";
		warn join ' ', @$c;
		warn "\n";
		warn "data points are:\n";
		warn join ' ', @$d;
		warn "\n";
		die "Error: columns don't match up!\n";
	}

	my %data;
	for my $idx ( 0 .. scalar(@$d) - 1 ) {
		$data{ $c->[$idx] } = $d->[$idx];
	}

	return \%data;
}

sub trim_whitespace {
	my $str = shift;
	my $copy = $str;
	chomp $copy;
	$copy =~ s/^\s+//g;
	$copy =~ s/\s+$//g;

	return $copy;
}

sub read_csts {
	my $fn = shift;
	open FILE, "<$fn" or die "Error opening file $fn ($!)";

	my %csts;
	while ( my $line = <FILE> ) {
		if ( $line =~ /SOGFUNC/ ) {
			$line = trim_whitespace( $line );
			my @d = split /\s+/, $line;

			my $resi = $d[2];
			my $resj = $d[4];
			for ( 1 .. 6 ) {
				shift @d;
			}

			shift @d;
			my @csts;
			while ( @d ) {
				my $mean   = shift @d;
				my $sdev   = shift @d;
				my $weight = shift @d;
				if ( $sdev == 0 ) {
					$sdev = 2.0;
				}
				my $cst = { mean => $mean, sdev => $sdev, weight => $weight };
				push @{ $csts{$resi}{$resj} }, $cst;
				#use Data::Dumper::Simple;
				#print "line = $line\n";
				#print Dumper($cst);
				#exit 1;
			}
			#if ( $resi == 49 && $resj == 65 ) {
			#	use Data::Dumper::Simple;
			#	print "before\n";
			#	print Dumper( $csts{$resi}{$resj} );
			#	@{ $csts{$resi}{$resj} } = @{ reweight_csts( $csts{$resi}{$resj} ) };
			#	print "after\n";
			#	print Dumper( $csts{$resi}{$resj} );
			#	exit 1;
			#}
			@{ $csts{$resi}{$resj} } = @{ reweight_csts( $csts{$resi}{$resj} ) };
		}
	}

	return \%csts;
}

sub merge_redundant_csts {
   my $csts = shift;

   my %merged_csts;
   foreach my $cst (@$csts) {
      $merged_csts{$cst->{mean}}{$cst->{sdev}} += $cst->{weight};
   }

	use constant min_cst_weight => 1e-5;
   my @new_csts;
   foreach my $mean ( keys %merged_csts ) {
   foreach my $sdev ( keys %{ $merged_csts{$mean} } ) {
      my $weight = $merged_csts{$mean}{$sdev};
		if ( $weight > min_cst_weight ) {
      	push @new_csts, { mean => $mean, sdev => $sdev, weight => $weight };
		}
   }
   }

   return \@new_csts;
}

sub reweight_csts {
	my $csts = shift;

	use constant k => 10;
	use constant min_cst_weight => 1e-5;
	my $total = 0;
	my $bg_weight = 0;
	foreach my $cst (@$csts) {
		if ( $cst->{mean} > 12 ) {
			$bg_weight += $cst->{weight};
		} else {
			$total += $cst->{sdev} ** (-1 * k)
		}
	}

	my $total_cst_count = scalar(@$csts);
	my $inf_wt = 1;
	if ( $bg_weight > 0 ) {
		$inf_wt = 1 - $bg_weight;
	}

	my @new_csts;
	foreach my $cst (@$csts) {
		my $new_cst = $cst;
		if ( $new_cst->{mean} <= 12 ) {
			$new_cst->{weight} = $inf_wt * $cst->{sdev} ** ( -1 * k ) / $total;
		}
		if ( $new_cst->{weight} > min_cst_weight ) {
			push @new_csts, $new_cst;
		}
		#push @new_csts, $new_cst;
	}
	#use List::Util qw/ sum /;
	#my $total_wt = sum map { $_->{weight} } @new_csts;
	#my $old_total = sum map { $_->{weight} } @$csts;

	#if ( $total_wt - 1 > 1e-5 ) {
	#	use Data::Dumper::Simple;

	#	print "total_wt  = $total_wt\n";
	#	print "old_total = $old_total\n";
	#	print "inf_wt = $inf_wt\n";

	#	foreach my $cst (@$csts) {
	#		my $new_cst = $cst;
	#		print "weight goes from ", $new_cst->{weight};
	#		if ( $new_cst->{mean} < 12 ) {
	#			$sdev = $cst->{sdev};
	#			print " ( $inf_wt * ", $cst->{sdev}, " ** ( -1 * 10 ) / $total ) ";
	#			$new_cst->{weight} = $inf_wt * $cst->{sdev} ** ( -1 * k ) / $total;
	#		}
	#		print " to ", $new_cst->{weight}, "\n";
	#	}
	#	print Dumper( $csts );
	#	print Dumper( @new_csts );
	#	die "Error: total_wt = $total_wt!\n";
	#}

	return \@new_csts;
}

sub dgaussian {
	my $mean  = shift;
	my $sdev  = shift;
	my $value = shift;

	use List::Util qw/ max /;
	use constant TINY_FLOAT => 1e-200;

	my $r = abs( $mean - $value );

	use constant sqrt_2pi => 2.50662721600161;
	my $val = 1 / (sqrt_2pi * $sdev) * exp( -1 * $r * $r / ( 2 * $sdev * $sdev ) );

	return max( $val, TINY_FLOAT );
}
