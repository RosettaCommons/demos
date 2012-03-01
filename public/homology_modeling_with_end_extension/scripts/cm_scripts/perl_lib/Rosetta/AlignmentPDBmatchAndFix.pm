package Rosetta::AlignmentPDBmatchAndFix;

use File::Basename;

use Rosetta::Template;
use Rosetta::Job;
use Rosetta::Util;

sub make_alignmentMatchPDB {
	my $aln_fn       = shift;
	my $out_aln_fn   = shift;
	my $template_dir = shift;
	my $options      = shift;

	my $dir = $options->{aln_dir};
	my $aln_fn_local = assemble_path($dir,$aln_fn);
	my @template_ids = @{ Rosetta::Template::ids_from_aln_file($aln_fn_local)};
	my @template_pdbs;
	my %seen;
	foreach my $id (@template_ids) {
		my $pdbid = lc( substr( $id, 0, 4 ) );
		my $chain = uc( substr( $id, 4, 1 ) );
		my $pdb_file = "$pdbid$chain.pdb";
		my $pdb_location = assemble_path( $template_dir, $pdb_file );
		push(@template_pdbs,$pdb_location) unless $seen{$pdb_file}++;
	}

	my $instance = Rosetta::Instance->new(
		mini_prefix => $options->{mini_path},
		compiler    => $options->{mini_compiler},
		mode        => $options->{mini_compile_mode},
		db_path     => $options->{mini_db_path},
	);
	my $job = $instance->generate_job( "fix_alignment_to_match_pdb" );
	$job->add_args(
		"-in:file:alignment $aln_fn",
		"-in:file:template_pdb ", @template_pdbs,
		"-database ", $options->{mini_db_path},
		"-out:file:alignment $out_aln_fn",
		"-cm:aln_format grishin",
		       );
	$job->dir($dir);
	$job->lockfile( "$aln_fn.alignmentPDB.lock" );
	$job->results_file($out_aln_fn);
	return $job;
}

1;
