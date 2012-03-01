package Rosetta::Constraints;

use File::Basename;
use Rosetta::Util;

sub make_constraints {
	my $aln_fn      = shift;
	my $fasta_fn    = shift;
	my $ev_map_file = shift;
	my $cst_file    = shift;
	my $options     = shift;

	my $dir = $options->{aln_dir};

	my $script_dir = Rosetta::Util::script_dir();
	my $job = Rosetta::Job->new(
		executable => assemble_path( $script_dir, 'predict_distances.pl' ),
		args => [
			$aln_fn, basename($fasta_fn),
			"--outfile $cst_file",
			"--ev_map_file $ev_map_file",
			"--aln_format grishin",
			"--max_templates", $options->{max_templates},
			"--max_e_value", 9999999,
		],
		dir => $dir,
		lockfile => "$cst_file.make_csts.lock",
		results_file => $cst_file
	);

	if ( exists $options->{template_dir} && -d $options->{template_dir} ) {
		$job->add_args( '--pdb_dir ', $options->{template_dir} );
	}

	return $job;
}

1;
