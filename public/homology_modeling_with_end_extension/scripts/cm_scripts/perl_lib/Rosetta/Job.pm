package Rosetta::Job;

use Cwd;
use File::Path;
use Exporter;

@ISA = qw/ Exporter /;
@EXPORT_OK = qw/ run_jobs_in_parallel /;

sub new {
	my $class = shift;

	# always need values these fields for internal consistency.
	my @required_fields = qw/ func_ref executable ok_without_results /;
	my %data = map { $_ => 0 } @required_fields;
	my $self = \%data;

	bless $self, $class;
	while ( @_ ) {
		my $key = shift;
		my $val = shift;
		if ( $key eq 'args' ) {
			$self->args( $val );
		} else {
			$self->{$key} = $val;
		}
	}

	return $self;
}

sub lockfile {
	my $self = shift;
	if ( @_ ) {
		$self->{lockfile} = shift;
	}
	return $self->{lockfile};
}

sub mpi {
	my $self = shift;
	if ( @_ ) {
		$self->{mpi} = shift;
	}
	return $self->{mpi};
}

sub results_file {
	my $self = shift;
	if ( @_ ) {
		$self->{results_file} = shift;
	}
	return $self->{results_file};
}

sub output_file {
	my $self = shift;
	if ( @_ ) {
		$self->{output_file} = shift;
	}
	return $self->{output_file};
}

sub executable {
	my $self = shift;
	if ( @_ ) {
		$self->{executable} = shift;
	}
	return $self->{executable};
}

sub func_ref {
	my $self = shift;

	if ( @_ ) {
		$self->{func_ref} = shift;
	}
	return $self->{func_ref};
}


sub args {
	my $self = shift;
	if ( @_ ) {
		$self->{args} = shift;
	}

	# Internally we want args to be stored as an arrayref, but sometimes people
	# want to pass in a single string.
	if ( ref( $self->{args} ) ne 'ARRAY' ) {
		$self->{args} = [ $self->{args} ];
	}

	return $self->{args};
}

sub add_args {
	my $self = shift;

	while ( @_ ) {
		my $arg = shift;
		push @{ $self->{args} }, $arg;
	}

	return $self;
}

sub clear_args {
	my $self = shift;
	$self->{args} = [];
}

sub args_as_str {
	my $self = shift;
	my $str  = join ' ', @{ $self->args };
	return $str;
}

sub output {
	my $self = shift;
	if ( @_ ) {
		$self->{output} = shift;
	}

	return $self->{output};
}

sub logfile {
	my $self = shift;

	if ( @_ ) {
		$self->{logfile} = shift;
	}

	# hack to always log output!
	if ( !exists $self->{logfile} ) {
		return join '.', ( $self->results_file, 'log' );
	}
	return $self->{logfile};
}

sub dir {
	my $self = shift;
	if ( @_ ) {
		$self->{dir} = shift;
	}
	return $self->{dir};
}

sub overwrite {
	my $self = shift;
	if ( @_ ) {
		$self->{overwrite} = shift;
	}

	return $self->{overwrite};
}

sub make_flags_file {
	my $self = shift;
	my $fn   = shift;

	if ( !$fn ) {
		die "Error: need a filename in method make_flags_file!\n";
	}
	if ( -f $fn ) {
		die "Error: not overwriting fn $fn!\n";
	}

	open FILE, ">$fn" or die $!;
	my @args = @{ $self->args };
	print FILE join "\n", @args;
	print FILE "\n";
	close FILE or die $!;
}

sub to_string {
	my $self = shift;
	my $str = "JOB:\n";

	my %data;
	foreach my $key ( keys %{ $self } ) {
		if ( ref($self->{$key}) eq 'ARRAY' ) {
			$data{$key} = join ' ', @{ $self->{$key} };
		} else {
			$data{$key} = $self->{$key};
		}
	}

	$str .= join "\n",
		map { "\t$_: " . $data{$_} } keys %data;

	return $str;
}

sub from_string {
	my $self = shift;
	my $str  = shift;

	my @lines = split /\n/, $str;
	foreach my $line (@lines) {
		if ( $line !~ /JOB/ && $line !~ /^\s*$/ ) {
			$line =~ s/^\s+//g;
			$line =~ s/\s+$//g;

			my ($key,@data) = split /\s+/, $line;
			$key =~ s/://g;
			if ( scalar(@data) == 1 ) {
				$self->{$key} = $data[0];
			} else {
				$self->{$key} = \@data;
			}
		}
	}
}

sub run_with_message {
	my $self = shift;
	my $message = shift;

	print "$message ... ";
	if ( $self->results_completed ) {
		print "already finished.\n";
	} elsif ( !$self->can_run ) {
		print "can't run.\n";
	} else {
		$self->run;
		print "done.\n";
	}
	return;
}

sub run {
	my $self = shift;

	if ( !$self->executable && !$self->func_ref ) {
		die "Error: no executable or function reference defined!\n";
	}
	if ( !$self->args ) {
		die "Error: no args defined!\n";
	}
	my $run = 1;
	my $lockfile = $self->lockfile;
	my $results_file = $self->results_file;

	if ( $self->overwrite ) {
		#print "Redoing job ...\n";
		#print "Removing lockfile $lockfile and result $results_file!\n";
		unlink_or_print( $lockfile );
		unlink_or_print( $results_file );
	}

	if ( $self->is_locked ) {
		print "Skipping run because job is locked with $lockfile!\n";
		#return;
		$run = 0;
	}

	if ( $self->results_file && $self->results_completed ) {
		print "Skipping run because $results_file looks completed!\n";
		$run = 0;
		#return;
	}

	my $original_dir = getcwd;
	if ( $self->dir && ! -d $self->dir ) {
		mkpath( $self->dir );
	}

	if ( $self->dir ) {
		#$cmd = "cd " . $self->dir . "; $cmd";
		#print "changing to directory ", $self->dir, "\n";
		chdir( $self->dir );
	}

	if ( $run ) {
		if ( !$self->acquire_lock ) {
			#print "changing to directory ", $original_dir, "\n";
			chdir $original_dir;
			return;
		}

		if ( $self->func_ref ) {
			# execute function reference, be sure to dereference
			# arguments into an array.
			my $func_ref = $self->func_ref;

			my @args = @{ $self->args };
			$func_ref->( @{ $self->args } );
		} elsif ( $self->executable ) {
			my $cmd = join ' ', ( $self->executable, $self->args_as_str );

			if ( $self->{mpi} ) {
				$cmd =~ s/\.default\.linuxgccrelease/\.mpi\.linuxgccrelease/g;
				$cmd =~ s/\.linuxgccrelease/\.mpi\.linuxgccrelease/g;
				$cmd = join ' ', ( $cmd, '-jd2:mpi_filebuf_jobdistributor' );
				$cmd = "smart_mpi.pl -cmd \'$cmd\'";
			}

			#print "executing cmd ", $cmd, "\n";
			my $output = `$cmd 2>&1`;
			$self->output( $output );
			if ( $self->logfile ) {
				my $logfile = $self->logfile;
				if ( -f $logfile && !$self->overwrite ) {
					print "Not overwriting output into $logfile!\n";
				} else {
					use File::Path;
					use File::Basename;
					mkpath( dirname($logfile) );
					open  FILE, ">$logfile" or die "Error opening $logfile ($!)\n";
					print FILE $output, "\n";
					close FILE or die $!;
				}
			}
		}
		$self->release_lock;
	}

	
	if ( $run && ! $self->results_completed ) {
		if ( !$self->{ok_without_results} ) {
			print "currently in ", getcwd, "\n";
			die "Error: no results from job!\n" . $self->to_string . "\n";
		}
	}

	if ( $self->dir ) {
		#print "changing to directory ", $original_dir, "\n";
		chdir $original_dir;
	}

}

sub can_run {
	my $self = shift;
	if ( $self->overwrite ) {
		return 1;
	}

	#if ( $self->is_locked ) {
	#	print "locked!\n";
	#}
	#if ( $self->results_completed ) {
	#	print "results_completed!\n";
	#}

	return ( !($self->is_locked || $self->results_completed) );
	return 1;
}

sub unlink_or_print {
	my $to_unlink = shift;

	my @files;
	if ( ref $to_unlink eq 'ARRAY' ) {
		@files = @{ $to_unlink };
	} elsif ( ref $fn eq 'HASH' ) {
		@files = @{ keys %$to_unlink };
	} else {
		push @files, $to_unlink;
	}

	foreach my $fn (@files) {
		if ( -f $fn ) {
			unlink $fn or print "Error unlinking file $fn!\n";
		}
	}
}

sub exit_if_not_finished {
	my $self = shift;
	if ( ! $self->finished ) {
		print "Exiting because job isn't finished!\n";
		print $self->to_string, "\n";
		exit 1;
	}
}

sub finished {
	my $self = shift;
	return ( !$self->is_locked && $self->results_completed );
}

sub results_completed {
	my $self = shift;

	if ( ! defined $self->results_file || ! $self->results_file ) {
		return 1;
	}

	my $dir = $self->dir;
	#if ( ! -d $dir ) {
	#	return 0;
	#}
	#my $orig_dir = getcwd;

	#chdir( $dir );
	my $results_file = $self->results_file;
	my $retval = 0;
	if ( $results_file && ( -f $results_file || -d $results_file || -f "$dir/$results_file" ) ) {
		$retval = 1;
	}
	#chdir( $orig_dir );

	return $retval;
}

sub is_locked {
	my $self = shift;

	my $lockfile = $self->lockfile;
	if ( !$lockfile ) {
		die "Error: no lockfile defined in job!\n";
	}
	return ( -f $lockfile );
}

sub acquire_lock {
	my $self = shift;

	if ( $self->is_locked ) {
		#print "Error: couldn't lock!\n";
		#print $self->to_string, "\n";
		return 0;
	}

	my $lockfile = $self->lockfile;
	system( "touch $lockfile" );
	return 1;
}

sub release_lock {
	my $self = shift;

	my $lockfile = $self->lockfile;
	if ( ! $self->is_locked ) {
		print "Error: trying to release lock but $lockfile doesn't exist!";
		return;
	}

	unlink( $lockfile );
}

sub print_logger {
	for my $msg (@_) {
		print $msg, "\n";
	}
}

sub run_jobs_in_parallel {
	my $jobs    = shift;
	my $n_procs = shift;

	my $logger = \&print_logger;
	if ( @_ ) {
		$logger = shift;
	}

	use Parallel::ForkManager;
	my $pm = Parallel::ForkManager->new($n_procs);
	for my $job (@$jobs) {
		$pm->start and next;
		my $PID = $$;
		my $results_file = $job->results_file;
		if ( $job->can_run ) {
			$logger->( "starting process $PID to make $results_file." );
			$job->run;
			$logger->( "finished process $PID." );
		}
		$pm->finish;
	}
	$pm->wait_all_children;
}

1;
