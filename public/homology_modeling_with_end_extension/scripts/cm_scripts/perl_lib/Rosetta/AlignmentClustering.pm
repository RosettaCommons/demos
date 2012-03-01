package Rosetta::AlignmentClustering;

use File::Basename;

use Rosetta::Template;
use Rosetta::Job;
use Rosetta::Util;

sub make_alignmentClusters {
	my $aln_fn       = shift;
	my $fasta_fn     = shift;
	my $template_dir = shift;
	my $ev_map_file  = shift;
	my $hh_map_file  = shift;
	my $options      = shift;

	my $dir = $options->{aln_dir};
	my $aln_fn_local = assemble_path($dir,$aln_fn);
	$fasta_fn = basename $fasta_fn;
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
	my $results_file = "alignmentCluster_1.filt"; # should make at least one cluster

	my $instance = Rosetta::Instance->new(
		mini_prefix => $options->{mini_path},
		compiler    => $options->{mini_compiler},
		mode        => $options->{mini_compile_mode},
		db_path     => $options->{mini_db_path},
	);
	my $job = $instance->generate_job( "cluster_alns" );
	$job->add_args(
		"-in:file:alignment $aln_fn",
		"-in:file:fasta $fasta_fn",
		"-cm:aln_format grishin",
		"-ev_map ev_map.txt",
		"-hh_map hh_map.txt",
		"-in:file:template_pdb ", @template_pdbs,
		"-database ", $options->{mini_db_path},
	);
	if ( `hostname` !~ /nrb/ ) {
		$job->add_args("-run:debug"),
	}
	$job->dir($dir);
	$job->lockfile( "$aln_fn.make_clstAlns.lock" );
	$job->results_file( $results_file );
	return $job;
}

1;
