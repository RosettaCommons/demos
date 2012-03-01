package Rosetta::Sequence;

use Exporter;
@ISA = qw/ Exporter /;
@EXPORT_OK = qw/
	pssm_from_fasta
	read_protein_fasta_fn
	write_fasta_file
        filter_alns_by_date
        filter_alns_by_date_job
/;

use Cwd;
use File::Basename;

use alignment;
use Rosetta::AA;
use Rosetta::Util;
use Rosetta::Sequence;

sub pssm_from_fasta {
	my $fasta = shift;
	my $dir   = shift;

	mkdir_safe( $dir );
	copy_safe( $fasta, $dir );
	my $fn = basename $fasta;
	my $profile_job = Rosetta::Job->new(
		dir          => $dir,
		args         => $fn,
		lockfile     => "$fn.profile.lock",
		executable   => "/work/tex/bin/run-psiblast.pl",
		results_file => "$fn.1.psiblast",
	);

	return $profile_job;
}

sub filter_alns_by_date_job {
	my $aln_fn            = shift;
	my $filtered_aln_file = shift;
	my $dir               = shift;
	my $date_cutoff       = shift;
	my $options           = shift;

	my $job = Rosetta::Job->new(
		func_ref     => \&filter_alns_by_date,
		args         => [ $aln_fn, $filtered_aln_file, $dir, $date_cutoff, $options ],
		dir          => $dir,
		date_cutoff  => $date_cutoff,
		lockfile     => "filter.lock",
		results_file => $filtered_aln_file,
	);

	return $job;
}

sub filter_alns_by_date {
	my $aln_fn            = shift;
	my $filtered_aln_file = shift;
	my $dir               = shift;
	my $date_cutoff       = shift;
	my $options           = shift;

	use Rosetta::PDB_Date;

	my $pdb_date = Rosetta::PDB_Date->new();

	if ( -f assemble_path( $dir, $filtered_aln_file ) ) {
		print "Skipping run because $filtered_aln_file exists!\n";
		return;
	}

	my $start_dir = getcwd;
	chdir( $dir );
	my @alns = @{ alignment::parse_alignments( $aln_fn ) };

	my %rejected;
	my $n_printed = 0;
	open FILE, ">$filtered_aln_file" or die
		"Error opening file $filtered_aln_file ($!)";

	ALN: foreach my $aln (@alns) {
		my $id = lc substr( $aln->get_template_name, 0, 4 );
		my $date = $pdb_date->date( $id );

		# if $date_cutoff is >0, then use it to try and
		# filter alignments by date.
		my $accept = 1;
		if ( $date_cutoff ) {
		    if ( $date && $date < $date_cutoff ) {
				$accept = 1;
			} else {
				$accept = 0;
			}
		}

		if ( $accept ) {
			$aln->template_name(
				lc( substr( $aln->get_template_name, 0, 4 ) ) .
				uc( substr( $aln->get_template_name, 4, 1 ) ) .
				  ( substr( $aln->get_template_name, 5 )    )
			);
			print FILE $aln->filt_string;
			$n_printed++;
		} else {
			$rejected{$id}++;
			print "rejected alignment to $id";
			print " (date = $date, cutoff = $date_cutoff)\n";
		}

		#if ( $n_printed >= $options->{max_templates} ) {
		#	print "stopped looking for alignments\n";
		#	print "picked $n_printed, max_templates = ",
		#		$options->{max_templates}, "\n";
		#	last ALN;
		#}
	}
	close FILE or die $!;

	foreach my $id ( sort keys %rejected ) {
		my $count = $rejected{$id};
		my $date = $pdb_date->date( $id );
		print join ' ', (
			"rejected $count alignments to", $id, "with date",
			$date, "(cutoff = $date_cutoff)"
		);
		print "\n";
	}

	print "$filtered_aln_file: selected ", $n_printed, " alignments from ",
		scalar(@alns), "\n";

	chdir( $start_dir );
}

sub remove_redundant_alns {
	my $aln_in  = shift;
	my $aln_out = shift;
	my $options = shift;

	use strict;

	my $file_in  = assemble_path( $options->{aln_dir}, $aln_in );
	my $file_out = assemble_path( $options->{aln_dir}, $aln_out );

	my @alns = @{ alignment::parse_alignments($file_in) };
	my @alns_out;

	print "examining ", scalar(@alns), " alignments.\n";
	if ( $options->{max_template_pct_id} >= 1 ) {
		@alns_out = @alns;
	} else {
		@alns_out =
			grep { $_->get_percent_identity <= $options->{max_template_pct_id} }
			@alns;
	}

	open FILE, ">$file_out" or die $!;
	foreach my $aln (@alns_out) {
		print FILE $aln->filt_string, "\n";
	}
	close FILE or die $!;
}

sub make_hhsearch_e_value_map {
	my $path = shift;
	my $out_fn = shift;
	my $max_e_value = shift;

	if ( -f assemble_path($path,$out_fn) && !-z assemble_path($path,$out_fn) ) {
		return;
	}

	use hhsearch;
	use List::Util qw/ min /;

	my @alns;
	my @files = glob( "$path/hhsearch_result_*" );
	push @files, glob( "$path/*.hhr" );
	push @alns, map { @{ hhsearch::parse_alignments($_) } } @files;

	if ( scalar(@files) == 0 ) {
		# hack - try parsing e-values from BLAST results
		# if hhsearch results don't exist
		@files = glob( "$path/*.blast" );
		push @alns, map { @{ blast::parse_alignments($_) } } @files;
	}

	if ( -f "$path/$out_fn" ) {
		return;
	}

	my %data;
	my @errors;
	#foreach my $fn (@files) {
		#my @alns = @{ hhsearch::parse_alignments( $fn ) };
		#print "read ", scalar(@alns), " from $fn\n";
	foreach my $aln (@alns) {
		my $template_name = $aln->template_name;
		my $e_value = $aln->get_e_value;

		if ( exists $data{$template_name} ) {
			my $err = abs( $data{$template_name} - $e_value );
			#print join ' ', ( $template_name, $err );
			#print "\n";
			push @errors, $err;
		}


		if ( $e_value < $max_e_value ) {
			if ( exists $data{$template_name} ) {
				if ( $e_value < $data{$template_name} ) {
					$data{$template_name} = $e_value;
				}
			} else {
				$data{$template_name} = $e_value;
			}
		}
	}

	open FILE, ">$path/$out_fn" or die $!;
	print FILE join ' ', qw/ template e_value /;
	print FILE "\n";
	foreach my $templ ( keys %data ) {
		my $e_value = $data{$templ};
		print FILE join ' ', ( $templ, $e_value );
		print FILE "\n";
	}
	close FILE or die $!;
}

sub make_hhsearch_e_value_map_job {
	my $path = shift;
	my $out_fn = shift;
	my $max_e_value = shift;

	my $job = Rosetta::Job->new(
		args         => [ $path, $out_fn, $max_e_value ],
		lockfile     => 'ev_map.lock',
		dir          => $path,
		results_file => $out_fn,
		func_ref     => \&make_hhsearch_e_value_map,
	);

	return $job;
}

sub write_fasta_file {
	my $id  = shift;
	my $seq = shift;
	my $fn  = shift;

	if ( -f $fn ) {
		die "Error: not writing to file $fn as it already exists!\n";
	}
	use File::Basename qw/ dirname /;
	if ( ! -d dirname($fn) ) {
		die "Error: directory for output file $fn doesn't exist!\n";
	}

	$seq =~ s/^>.*\n//g; # remove old id

	open FILE, ">$fn" or die $!;
	print FILE "> $id\n";
	print FILE "$seq\n";
	close FILE or die $!;
}

sub read_fasta_fn {
	my $fasta_fn = shift;
	open FILE, "<$fasta_fn" or die "Error opening file $fasta_fn ($!)";
	my @file = <FILE>;
	close FILE or die $!;

	return join '',
		map { chomp $_; $_ }
		grep { !/^>/ }
		map { $_ =~ s/\s+//g; $_ }
		@file;
}

sub read_protein_fasta_fn {
	my $fasta_fn = shift;
	my $seq = read_fasta_fn($fasta_fn);

	# do a coherency check on the fasta sequence to make sure that all letters
	# can plausibly represent amino acids.
	my $idx = 1;
	foreach my $aa (split //, $seq) {
		if ( !Rosetta::AA::is_valid_name1($aa) ) {
			print "Error with sequence from $fasta_fn!\n";
			print "\t", "$aa at position $idx doesn't look like an amino acid!\n";
		}
		$idx++;
	}

	return $seq;
}

sub assign_unique_names {
	my $alns = shift;

	my %names;
	foreach my $aln (@$alns) {
		my $count = 1;
		my @d = split /_/, $aln->template_name;
		$aln->template_name( join '_', ( $d[0], $count ) );

		while ( exists $names{$aln->template_name} ) {
			$count++;
			$aln->template_name( join '_', ( $d[0], $count ) );
		}
		$names{$aln->template_name}++;
	}
}

sub convert_aln {
	my $format_in  = shift;
	my $format_out = shift;
	my $file_in    = shift;
	my $file_out   = shift;
	my $options    = shift;

	use strict;

	my $max_pct_id = 1.0;
	if ( exists $options->{max_pct_id} ) {
		$max_pct_id = $options->{max_pct_id};
	}

	print "converting $file_in ($format_in) => $format_out ($format_out)\n";

	my @alns;
	if ( $format_in eq 'blast' ) {
		@alns = @{ blast::parse_alignments($file_in) };
	} elsif ( $format_in eq 'grishin' ) {
		@alns = @{ alignment::parse_alignments($file_in) };
	} elsif ( $format_in eq 'hhsearch' ) {
		@alns = @{ hhsearch::parse_alignments($file_in) };
	}

	# give all of the alignments unique names.
	assign_unique_names(\@alns);

	print "read ", scalar(@alns), " alignments from $file_in.\n";
	open FILE, ">$file_out" or die $!;
	foreach my $aln (@alns) {
		my $pct_id = sprintf( "%2.2f", $aln->get_percent_identity );
		if ( $pct_id <= $max_pct_id ) {
			$aln->template_name(
				lc( substr( $aln->get_template_name, 0, 4 ) ) .
				uc( substr( $aln->get_template_name, 4, 1 ) ) .
					( substr( $aln->get_template_name, 5 )   )
      	);
      	#my $id = substr( $aln->template_name, 0, 5 );
			print "accepted alignment with id ", $aln->template_name;
			print " ($pct_id <= $max_pct_id)\n";

			print FILE $aln->filt_string;
		}
	}
	close FILE or die $!;
}

1;
