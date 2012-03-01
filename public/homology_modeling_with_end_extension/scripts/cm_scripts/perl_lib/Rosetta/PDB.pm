package Rosetta::PDB;
require Exporter;

@ISA = qw/ Exporter /;
@EXPORT_OK = qw/
	pdb_to_seq
	pdb_to_fasta
/;

use alignment;

use Rosetta::Job;
use Rosetta::AA;
use Rosetta::Sequence;
use Rosetta::Util;
use Rosetta::Instance;

use File::Basename;

sub idealize_pdb {
	my $fn = shift;
	my $dest_dir = shift;
	my $options  = shift;

	copy_safe( $fn, $dest_dir );
	my $pdb = basename $fn;

	my $results_file = $pdb;
	$results_file =~ s/\.pdb$//g;
	$results_file = join '', ($results_file,'_0001','.pdb');

	my $instance = Rosetta::Instance->new(
		mini_prefix => $options->{mini_path},
		mode        => $options->{mini_compile_mode},
		compiler    => $options->{mini_compiler}
	);

	my $job = $instance->generate_job(
		"idealize",
		"-in::file::s $pdb -database $options->{mini_db_path}",
	);
	$job->dir( $dest_dir );
	$job->lockfile( "$pdb.idealize.lock" );
	$job->results_file( $results_file );

	return $idealize_job;
}

sub pdb_to_fasta {
	my $pdb_fn = shift;

	my $base = basename $pdb_fn;
	my $seq = "> $base sequence\n";
	$seq .= pdb_to_seq($pdb_fn);

	return "$seq\n";
}

sub pdb_to_seq {
	my $pdb_fn = shift;

	#print "opening ($pdb_fn)\n";
	open FILE, "<$pdb_fn" or die "Error opening file $pdb_fn ($!)";
	my $last_resnum = 0;
	my $seq = '';
	while ( my $line = <FILE> ) {
		if ( length($line) > 20 && $line =~ /^ATOM/ ) {
			my $resnum = substr( $line, 22, 4 );
			if ( $resnum ne $last_resnum ) {
				my $name3 = substr( $line, 17, 3 );
				$seq .= Rosetta::AA::name3_to_name1( $name3 );
			}
			$last_resnum = $resnum;
		}
	}

	close FILE or die $!;
	return $seq;
}

1;
