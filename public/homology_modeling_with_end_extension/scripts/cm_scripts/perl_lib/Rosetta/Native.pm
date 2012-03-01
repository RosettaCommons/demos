package Rosetta::Native;

use Cwd qw/ abs_path /;
use Rosetta::Job;
use Rosetta::Util;
use Rosetta::Instance;
use File::Basename;

sub setup_native_job {
	my $options    = shift;
	my $native_dir = shift;

	my $final_native_fn = assemble_path( $native_dir, 'native.pdb' );

	my $job = Rosetta::Job->new(
		func_ref     => \&setup_native,
		args         => [ $options, $native_dir ],
		lockfile     => "setup_native.lock",
		dir          => $native_dir,
		results_file => $file_native_fn,
	);

	return $job;
}

sub full_length_job {
	my $options = shift;
	my $native_dir = shift;

	my $full_fasta = abs_path($options->{fasta_fn});
	my $native_fn  = abs_path($options->{native_fn});
	my $frag3_fn   = abs_path($options->{frag3_fn});

	my @args = (
		"-in:file:fasta $full_fasta",
		"-cm:aln_format grishin",
		"-cm:min_loop_size 5",
		"-loops:frag_files $frag3_fn none",
		"-loops:frag_sizes 3 1",
		"-in:file:s $native_fn",
	);

	my $instance = Rosetta::Instance->new(
		mini_prefix => $options->{mini_path},
		mode        => $options->{mini_compile_mode},
		compiler    => $options->{mini_compiler},
		db_path     => $options->{mini_db_path},
	);
	my $output_fn = join '_', ( $native_fn, "full_length.pdb" );
	my $full_length_job = $instance->generate_job( "full_length_model", join ' ', @args );
	$full_length_job->results_file( $output_fn );
	$full_length_job->lockfile( "full_length.lock" );
	$full_length_job->dir( $native_dir );

	return $full_length_job;
}

sub setup_native {
	my $options    = shift;
	my $native_dir = shift;

	my $full_fasta = $options->{fasta_fn};
	my $native_fn  = $options->{native};
	my $frag3_fn   = $options->{frag3_fn};

	if ( !$native_fn || ! -f $native_fn ) {
		print "No native given!\n";
		return 0;
	}
	#print "creating full_length native for $native_fn in $native_dir\n";

	mkdir_safe( $native_dir );
	copy_safe( $full_fasta, $native_dir );
	copy_safe( $native_fn, $native_dir );
	$native_fn = basename $native_fn;

	my @args = (
		"-in:file:fasta $full_fasta",
		"-cm:aln_format grishin",
		"-cm:min_loop_size 5",
		"-loops:frag_files $frag3_fn none",
		"-loops:frag_sizes 3 1",
		"-in:file:s $native_fn",
	);

	my $instance = Rosetta::Instance->new(
		mini_prefix => $options->{mini_path},
		mode        => $options->{mini_compile_mode},
		compiler    => $options->{mini_compiler},
		db_path     => $options->{mini_db_path},
	);
	my $output_fn = join '_', ( $native_fn, "full_length.pdb" );
	my $full_length_job = $instance->generate_job( "full_length_model", join ' ', @args );
	$full_length_job->results_file( $output_fn );
	$full_length_job->lockfile( "full_length.lock" );
	$full_length_job->dir( $native_dir );

	if ( ! -f $frag3_fn ) {
		die "Error: can't find fragment file $frag3_fn!\n";
	}

	$full_length_job->run_with_message( "making full-length native from $native_fn" );

	$native_fn = basename $native_fn;
	my $full_length_native = assemble_path( $native_dir,"full_$native_fn" );
	system( join ' ', ( "cp", assemble_path($native_dir,$output_fn), $full_length_native ) );

	$options->{native_fn} = assemble_path( $native_dir, "native.pdb" );
	symlink_safe(
		$full_length_native, $options->{native_fn}
	);
}

1;
