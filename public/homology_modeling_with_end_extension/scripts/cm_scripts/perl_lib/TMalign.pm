package TMalign;

use Exporter;

use constant tmalign_executable => '~tex/bin/TMalign';

@ISA = qw/ Exporter /;
@EXPORT_OK = qw/
	get_tmscore
/;

sub get_tmscore {
	my $model  = shift;
	my $native = shift;

	assert_files_exist($model,$native);

	# check for nan in coordinates
	if ( pdb_has_nans($model) ) {
		return 0.0;
	}

	my $cmd = join ' ', ( tmalign_executable, $model, $native );
	my $output = `$cmd`;

	if ( $output !~ /TM-score=([\d\.]+),/ ) {
		die "Error: can't find TM-score for ($model,$native) in output:\n$output\n";
	}
	my $tmscore = $1;
	return $tmscore;
}

sub assert_files_exist {
	foreach my $fn (@_) {
		if ( ! -f $fn ) {
			die "Error: file $fn doesn't exist!\n";
		}
	}
}

sub pdb_has_nans {
	my $fn = shift;

	my $has_nans = 0;
	open PDB, "<$fn" or die $!;
	while ( my $line = <PDB> ){
		if ( $line =~ /nan/ ) {
			$has_nans = 1;
		}
	}

	close PDB or die $!;
	return $has_nans;
}
