package Rosetta;

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
	);

	foreach my $key ( keys %defaults ) {
		if ( ! exists $self->{$key} ) {
			$self->{$key} = $defaults{$key};
		}
	}

	return $self;
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

sub run_program {
	my $self   = shift;
	my $binary = shift;
	my $args   = shift;

	my $mode     = $self->mode;
	my $compiler = $self->compiler;
	my $platform = $self->determine_platform;

	my $suffix     = join '' , ( $platform, $compiler, $mode );
	my $executable = join '.', ( $binary, $suffix );
	my $cmd        = join ' ', ( $executable, $args );
	my $output     = `$cmd`;

	return $output;
}

1;
