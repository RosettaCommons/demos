package Rosetta::TargetFiles;

use Rosetta::Util;

sub new {
	my $class       = shift;
	my $target_id   = shift;
	my $library_dir = shift;
	my $def_fn      = shift;

	if ( ! $target_id ) {
		die "Error: must pass (target_id,library_dir) for TargetFiles!\n";
	}

	if ( ! $library_dir || !-d $library_dir ) {
		die "Must define library_dir for TargetFiles ($library_dir)\n";
	}

	my $self = bless {}, $class;

	$self->value( 'target_id'   => $target_id );
	$self->value( 'library_dir' => $library_dir );
	$self->fill_in_defaults( $target_id, $library_dir, $def_fn );

	while ( @_ ) {
		my $key = shift;
		my $val = shift;
		$self->value( $key, $value );
	}

	return $self;
}

sub fill_in_defaults {
	my $self        = shift;
	my $target_id   = shift;
	my $library_dir = shift;
	my $def_fn      = shift;

	my %vars = (
		target_id   => $target_id,
		library_dir => $library_dir,
	);

	if ( !defined $def_fn || !-f $def_fn ) {
		print "Warning: file $def_fn doesn't exist!\n";
		$def_fn = assemble_path(
			cm_base(), "protocols", "protocol_files.conf",
		);
	}
	#print "Defining protocol files using $def_fn.\n";
	my %options;
	Rosetta::Util::load_options_from_file($def_fn,\%options);

	foreach my $key (keys %options) {
		$options{$key} = Rosetta::Util::eval_template_str(
			$options{$key}, \%vars
		);
		$self->value( $key, $options{$key} );
	}
}

sub value {
	my $self  = shift;
	my $key   = shift;

	if ( @_ ) {
		my $value = shift;
		$self->{$key} = $value;
		#print "setting $key to $value\n";
	}

	if ( !defined $self->{$key} ) {
		die "TargetFiles error: no option for $key\n" . $self->to_string;
	}

	return $self->{$key};
}

sub to_string {
	my $self = shift;

	my @keys = sort keys %{ $self };

	my $str = "def:\n";
	foreach my $key (@keys) {
		$str .= join " => ", ($key,$self->{$key});
		$str .= "\n";
	}

	return $str;
}

1;
