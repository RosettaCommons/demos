package Rosetta::Template;

use Exporter;
@ISA = qw/ Exporter /;

@EXPORT_OK = qw/
	ids_from_aln_file
	filter_alignments_by_id
/;

use Cwd;
use Rosetta::PDB;
use Rosetta::Job;
use Rosetta::Util;
use Rosetta::JobList;
use Rosetta::Sequence;

use blast;
use alignment;

use List::Util qw/ min /;
use File::Basename;

sub ids_from_aln_file {
	my $aln_fn = shift;
	my @alns = @{ alignment::parse_alignments( $aln_fn ) };

	my @template_ids = map { substr( $_, 0, 5 ) } map { $_->template_name }
		@alns;

	# make a unique list of template_ids
	my %count;
	my @unique;
	foreach my $id (@template_ids) {
		if ( !exists $count{$id} ) {
			$count{$id}++;
			push @unique, $id;
		}
	}

	return \@unique;
	#return unique( \@template_ids );
}

sub make_template_vall {
	my $ideal_pdb_path = shift;
	my $profile_path   = shift;

	my @ideal_pdbs = glob( "$ideal_pdb_path/*0001.pdb" );
	foreach my $pdb (@ideal_pdbs) {
		my $fn = basename $pdb;
		my $id = substr( $fn, 0, 5 );
		# read profile, ss and torsions
		my $profile = assemble_path( $profile_path, "$id.fasta.10.chk" );
	}
}

sub make_vall_lines {
	my $pdb_file   = shift;
	my $fasta_file = shift;
	my $dir        = shift;

	if ( ! -f $pdb_file ) {
		print STDERR "No PDB file found ($pdb_file)\n";
	}
	if ( ! -f $fasta ) {
		print STDERR "No fasta file found ($fasta)\n";
	}

	my $ideal_dir = assemble_path( $dir, 'ideal' );
	my $vall_dir  = assemble_path( $dir, 'vall' );
	my $pssm_dir  = assemble_path( $dir, 'pssm' );

	my $pssm_job = Rosetta::Sequence::pssm_from_fasta(
		assemble_path($pssm_dir,$fasta_file)
	);
	my $idealize_job = Rosetta::PDB::idealize_pdb(
		assemble_path($dir,$pdb_file), $ideal_dir
	);
	$pssm_job->run;
	$idealize_job->run;
	$vall_lines .= sprintf (
		"%5s %1s %1s %5d %4d %4d " .
		"%8.2f %8.2f %8.2f " .       # xyz
		"%8.3f %8.3f %8.3f %8.3f " . # torsions
		"%3d %4.2f %5.3f ",
		$pdb_id,
		$res_ids[$res_i],
		$ideal_info->[$i]->{ss},
		$real_res_i[$i],
		$protein_position_end,
		$protein_position_begin,
		$ideal_info->[$i]->{CA_x},
		$ideal_info->[$i]->{CA_y},
		$ideal_info->[$i]->{CA_z},
		$ideal_info->[$i]->{phi},
		$ideal_info->[$i]->{psi},
		$ideal_info->[$i]->{omg},
		$ideal_info->[$i]->{chi},
		$nalign,
		$acc,
		$gap
	);
}

sub make_template_alignments {
	my $fasta_fn = shift;
	my $aln_fn   = shift;
	my $dir      = shift;
	my $options  = shift;

	my $aln_maker = $options->{make_alignments};
	my $options_dump = '';
	if ( exists $options->{use_hhblits} ){
		$options_dump = "-use_hhblits";
	}
	mkdir_safe( $dir );
	copy_safe( $fasta_fn, $dir );
	my $job = Rosetta::Job->new(
		executable   => $aln_maker,
		args         => [
			basename($fasta_fn),
			'-outfile', $aln_fn,
			'-min_alns', 50, $options_dump,
			'-max_template_pct_id', $options->{max_template_pct_id},
			'-aln_dir', $options->{aln_dir},
			'-targetid', $options->{target_id},
		],
		lockfile     => "make_alns.lock",
		dir          => $dir,
		results_file => $aln_fn,
	);
	return $job;
}

sub setup_templates {
	my $template_dir    = shift;
	my $aln_file        = shift;
	my $options         = shift;
	my $native_fn       = shift;

	my $new_aln_fn = join '.', ( basename($aln_file), "valid" );
	if ( -f assemble_path( $options->{aln_dir}, $new_aln_fn ) ) {
		return;
	}

	$options->{logger}( "setting up templates in directory $template_dir" );

	mkdir_safe( $template_dir );

	my @valid_ids;
	my @template_ids = @{ ids_from_aln_file( $aln_file ) };
	my $n_templates  = min( $options->{max_templates}, scalar(@template_ids) );

	$options->{logger}( "validating templates!" );
	ID: foreach my $id (@template_ids) {
		my $pdbid = lc( substr( $id, 0, 4 ) );
		my $chain = uc( substr( $id, 4, 1 ) );
		my $pdb_file = "$pdbid$chain.pdb";
		#print "looking for $pdb_file\n";
		$options->{logger}( "looking for $pdb_file ... " );
		if (  -f assemble_path( $template_dir, $pdb_file ) ) {
			#print " already exists as ", assemble_path($template_dir,$pdb_file), "\n";
			$options->{logger}(
				" already exists as ", assemble_path($template_dir,$pdb_file)
			);
		} else {
			$options->{logger}( "attempting to fetch using get_pdb.py ... " );
			my $get_pdb = assemble_path( script_dir(), 'get_pdb.py' );

			my $cmd = "cd $template_dir; $get_pdb $pdbid $chain";
			my $output = `$cmd`;
			$options->{logger}($output);
			if ( ! -f assemble_path( $template_dir, $pdb_file ) ) {
				#print "Error getting pdb: $pdbid $chain\n";
				#print "output:\n$output\n";
				$options->{logger}( "Error getting pdb: $pdbid $chain" );
				$options->{logger}( "output:\n$output" );
			}
		}

		if ( ! -f assemble_path( $template_dir, $pdb_file ) ) {
			$options->{logger}( "error getting pdb: $pdbid $chain!" );
			next ID;
		} else {
			$options->{logger}( " found pdb for $pdbid $chain." );
		}

		my $valid_id = "$pdbid$chain";
		push @valid_ids, $valid_id;
		if ( scalar(@valid_ids) >= $n_templates ) {
			$options->{logger}( "stopping because I only want $n_templates alignments." );
			last ID;
		}
	} # foreach my $id (@template_ids)

	# only include alignments that have valid pdb files
	$options->{logger}( "filtering alignments with these ids:" );
	$options->{logger}( @valid_ids );
	filter_alignments_by_id(
		basename($aln_file), $new_aln_fn, dirname($aln_file), \@valid_ids
	);
	my $aln_dir = dirname($aln_file);

	unlink(assemble_path($aln_dir,'alignment.filt'));
	make_local_link( $aln_dir, $new_aln_fn, 'alignment.filt' );

}


sub filter_alignments_by_id {
	my $aln_fn = shift;
	my $filtered_aln_file = shift;
	my $dir = shift;
	my $ids_to_include = shift;

	if ( -f assemble_path( $dir, $filtered_aln_file ) ) {
		print "Skipping run because $filtered_aln_file exists!\n";
		return;
	}

	my $start_dir = getcwd;
	chdir( $dir );
	my @alns = @{ alignment::parse_alignments( $aln_fn ) };

	my %id_count;
	foreach my $id (@$ids_to_include) {
		$id_count{$id} = 0;
	}

	my $n_printed = 0;
	open FILE, ">$filtered_aln_file" or die
		"Error opening file $filtered_aln_file ($!)";
	foreach my $aln (@alns) {
		$aln->template_name(
			lc( substr( $aln->get_template_name, 0, 4 ) ) .
			uc( substr( $aln->get_template_name, 4, 1 ) ) .
			  ( substr( $aln->get_template_name, 5 )    )
		);
		my $id = substr( $aln->template_name, 0, 5 );

		if ( exists $id_count{$id} ) {
			print FILE $aln->filt_string;
			$id_count{$id}++;
		}
	}
	close FILE or die $!;

	foreach my $id ( sort keys %id_count ) {
		my $count = $id_count{$id};
		print join ' ', ( "included $count alignments to", $id );
		print "\n";
	}

	chdir( $start_dir );
}

1;
