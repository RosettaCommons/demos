package Rosetta::FileUtil;

@ISA = qw/ Exporter /;

@EXPORT_OK = qw/
	string_from_file
/;

sub string_from_file {
	my $fn = shift;

	my $str = '';
	open FILE, "<$fn" or die $!;
	while ( my $line = <FILE> ) {
		$str .= $line;
	}
	close FILE or die $!;

	return $str;
}

1;
