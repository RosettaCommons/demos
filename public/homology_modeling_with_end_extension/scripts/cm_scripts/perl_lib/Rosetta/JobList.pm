package Rosetta::JobList;

sub new {
	my $class = shift;

	my $self = {};
	while ( @_ ) {
		my $key = shift;
		my $val = shift;
		$self->{$key} = $val;
	}

	bless $self, $class;
	return $self;
}

sub jobs {
	my $self = shift;

	if ( @_ ) {
		$self->{jobs} = [];
		while ( @_ ) {
			$self->{jobs} = shift;
		}
	}

	return $self->{jobs};
}

sub run {
	my $self = shift;
	my $jobs = $self->{jobs};

	if ( !$self->can_run ) {
		return;
	}
	foreach my $job (@$jobs) {
		#print "running job:\n", $job->to_string, "\n";
		$job->run;
	}
}

sub n_jobs {
	my $self = shift;
	return scalar( @{ $self->jobs } );
}

sub results_file {
	my $self = shift;
	if ( @_ ) {
		$self->{results_file} = $_;
	}
	return $self->{results_file};
}

sub can_run {
	my $self = shift;

	if ( defined $self->results_file ) {
		#print "looking for ", $self->results_file, "\n";
		return ( !-f $self->results_file );
	}

	#foreach my $job (@{$self->jobs}) {
	#	print "checking if job can run:\n";
	#	use Data::Dumper::Simple;
	#	print Dumper( $job );
	#	print $job->to_string, "\n";
	#}

	my @problems =
		grep { !$_->can_run } @{ $self->jobs };

	#for my $job (@problems) {
	#	print "Can't run job:\n";
	#	print $job->to_string;
	#	print "\n";
	#}

	return( scalar(@problems) == 0 );
}

sub to_string {
	my $self = shift;

	foreach my $job (@{$self->jobs}) {
		print $job->to_string, "\n";
	}
}

sub run_with_message {
	my $self = shift;
	my $message = shift;

	print "$message ... ";
	if ( !$self->can_run ) {
		print "can't run.\n";
	}

	$self->run;
	print "done.\n";
	return;
}


1;
