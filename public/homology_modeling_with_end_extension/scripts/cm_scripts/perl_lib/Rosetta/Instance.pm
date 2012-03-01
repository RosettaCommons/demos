package Rosetta::Instance;

use Rosetta::Job;
use Rosetta::Util;

use constant DEFAULT_MINI_DB_PATH => $ENV{HOME} . '/minirosetta_database';

sub new {
	my $class = shift;

	my $self = {};
	while ( @_ ) {
		my $key = shift @_;
		my $val = shift @_;
		$self->{$key} = $val;
	}

	bless $self, $class;

	$self->set_defaults;
	return $self;
}

sub set_defaults {
	my $self = shift;

	my %defaults = (
		mode => 'release',
		compiler => 'gcc',
		mini_prefix => assemble_path( $ENV{HOME}, 'src/mini/bin' )
	);

	foreach my $key ( keys %defaults ) {
		if ( ! exists $self->{$key} ) {
			$self->{$key} = $defaults{$key};
		}
	}

	return $self;
}

sub mini_prefix {
	my $self = shift;

	if ( @_ ) {
		$self->{mini_prefix} = shift;
	}

	return $self->{mini_prefix};
}

sub determine_platform {
	my $uname_output = `uname`;

	my $platform = 'undefined';
	if ( $uname_output =~ /Darwin/ ) {
		$platform = 'macos';
	} elsif ( $uname_output =~ /Linux/ ) {
		$platform = 'linux';
	}

	return $platform;
}

sub compiler {
	my $self = shift;
	if ( @_ ) {
		$self->{compiler} = shift;
	}

	return $self->{compiler};
}

sub mode {
	my $self = shift;
	if ( @_ ) {
		$self->{mode} = shift;
	}

	return $self->{mode};
}

sub db_path {
	my $self = shift;
	if ( @_ ) {
		$self->{db_path} = shift;
	}

	# try to provide a reasonable database default if none exists
	# in the user's environment.
	if ( !exists $self->{db_path} ) {
		return DEFAULT_MINI_DB_PATH;
	}
	return $self->{db_path};
}

sub run_program {
	my $self   = shift;
	my $binary = shift;
	my $args   = shift;

	my $job        = $self->generate_job($binary,$args);
	my $executable = $job->executable;
	my $args       = $job->args;
	my $cmd        = join ' ', ( $executable, $args );
	my $output     = `$cmd`;

	return $output;
}

sub generate_job {
	my $self   = shift;
	my $binary = shift;
	my $args   = shift;

	my $mini_prefix = $self->mini_prefix;
	my $mode        = $self->mode;
	my $compiler    = $self->compiler;
	my $platform    = $self->determine_platform;

	# replace newlines with spaces
	$args =~ s/\n/ /g;
	if ( !exists $ENV{ROSETTA3_DB} ) {
		$args = join ' ', ( $args, "-database", $self->db_path );
	}

	my $suffix     = join '' , ( $platform, $compiler, $mode );
	my $executable = join '.', ( $binary, $suffix );

	if ( $mini_prefix ) {
		$executable = assemble_path( $mini_prefix, join '.', ( $binary, $suffix ) );
		if ( ! -f $executable ) {
			$executable = assemble_path( $mini_prefix, 'bin', join '.', ( $binary, $suffix ) );
		}
		if ( ! -f $executable ) {
			die "Error: executable $executable doesn't exist!\n";
		}
	}

	my $job = Rosetta::Job->new(
		executable   => $executable,
		args         => $args,
		dir          => $ENV{PWD},
	);
	return $job;
}

1;
