package Rosetta::Fragments;

use Rosetta::Job;
use Rosetta::JobList;

use Rosetta::Util;
use Rosetta::Sequence;

use File::Basename;

sub make_fragments {
	my $fasta    = shift;
	my $frag_dir = shift;
	my $options  = shift;

	my $target_id = $options->{target_id};

	my $temp_target_id  = 't001_';
	my $fasta_base = basename $fasta;
	if ( $fasta_base =~ /(.*)\.fasta/ && length($1) == 5 ) {
		$temp_target_id = $1;
	}

	mkdir_safe( $frag_dir );
	unlink( assemble_path( $frag_dir, $fasta ) );
	copy_safe( $fasta, $frag_dir );
	my $frag_job = Rosetta::Job->new(
		executable   => "~tex/src/nnmake/make_fragments.pl",
		#args         => "-noprof $fasta_base -id $temp_target_id",
		args         => "-nosam -nojufo -noprof $fasta_base -id $temp_target_id",
		dir          => $frag_dir,
		lockfile     => "frags.lock",
		logfile      => "make_fragments.log",
		results_file => join('', ("aa",$temp_target_id,"03_05.200_v1_3")),
	);

	my @fragment_files = (
		join('', ("aa",$temp_target_id,"03_05.200_v1_3")),
		join('', ("aa",$temp_target_id,"09_05.200_v1_3")),
	);

	my $rename_job = Rosetta::Job->new(
		func_ref     => \&symlink_files,
		args         => [ \@fragment_files, $temp_target_id, $target_id ],
		results_file => join('', ("aa",$target_id,"03_05.200_v1_3")),
		lockfile     => "rename_fragments.lock",
		dir          => $frag_dir,
	);
	my $joblist = Rosetta::JobList->new(
		name => "setup fragments",
		jobs => [ $frag_job, $rename_job ],
		results_file => assemble_path( $rename_job->dir, $rename_job->results_file ),
	);
	return $joblist;
}

sub symlink_files {
	my $files    = shift;
	my $orig_str = shift;
	my $new_str  = shift;

	foreach my $fn (@$files) {
		my $new_fn = $fn;
		$new_fn =~ s/$orig_str/$new_str/g;
		print STDERR "symlinking $fn -> $new_fn\n";
		#rename($fn, $new_fn);
		symlink_safe( $fn, $new_fn );
	}
}

sub boincify_fragments_in_path {
	my $path   = shift;
	my $prefix = shift;

	my @files = grep { basename($_) !~ /$prefix/ }
		glob( assemble_path( $path, '*v1_3' ) );
	#@files = glob( assemble_path( $path, '*v1_3' ) );
	foreach my $fn (@files) {
		my $dir  = dirname $fn;
		my $base = basename $fn;
		my $result_file = join '_', ($prefix,$base);
		my $job = Rosetta::Job->new(
			dir          => $dir,
			lockfile     => "$base.boincify.lock",
			results_file => $result_file,
			func_ref     => \&boincify_fragment_file,
			args         => [ $base, $result_file ],
		);
		if ( $job->can_run ) {
			$job->run_with_message( "boincifying $fn" );
		}

		#make_local_link( $dir, $result_file, basename $fn );
	}
}

sub boincify_fragment_file {
	my $fn        = shift;
	my $result_fn = shift;

	open FILE, "<$fn" or die "Error opening file $fn ($!)";
	open OUTFILE, ">$result_fn" or die "Error opening file $result_fn ($!)";

	while ( my $line = <FILE> ) {
		if ( length($line) > 50 ) {
			chomp $line;
			my $filtered_line = substr($line, 0, 44);
			print OUTFILE $filtered_line, "\n";
		} else {
			print OUTFILE $line;
		}
	}

	close FILE or die $!;
	close OUTFILE or die $!;
}

1;
