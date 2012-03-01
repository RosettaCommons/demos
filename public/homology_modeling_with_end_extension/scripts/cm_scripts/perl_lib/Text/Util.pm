package Text::Util;

@ISA = qw/ Exporter /;
@EXPORT_OK = qw/
	data_from_str
	read_data
	data_to_hashref
	trim_whitespace
/;

sub data_from_str {
	my $str = shift;

	my @lines = map { trim_whitespace($_) } split /\n/, $str;
	my $header = shift @lines;
	my @columns = split /\s+/, $header;

	my @data;
	foreach my $line (@lines) {
		my @d = split /\s+/, $line;
		push @data, data_to_hashref( \@d, \@columns );
	}
	return \@data;
}

sub read_data {
	my $fn = shift;

	my $cols_wanted = 0;
	if ( @_ ) {
		$cols_wanted = shift;
	}
	open FILE, "<$fn" or die $!;
	my $header = <FILE>;
	$header = trim_whitespace( $header );
	my @columns = split /\s+/, $header;
	my @data;
	while ( my $line = <FILE> ) {
		$line = trim_whitespace( $line );
		my @d = split /\s+/, $line;

		push @data, data_to_hashref( \@d, \@columns );
	}
	close FILE or die $!;

	return \@data;
}

sub data_to_hashref {
	my $d = shift;
	my $c = shift;

	if ( scalar(@$d) != scalar(@$c) ) {
		die "Error: columns don't match up!\n";
	}

	my %data;
	for my $idx ( 0 .. scalar(@$d) - 1 ) {
		$data{ $c->[$idx] } = $d->[$idx];
	}

	return \%data;
}

sub trim_whitespace {
	my $str = shift;
	my $copy = $str;
	chomp $copy;
	$copy =~ s/^\s+//g;
	$copy =~ s/\s+$//g;

	return $copy;
}

1;
