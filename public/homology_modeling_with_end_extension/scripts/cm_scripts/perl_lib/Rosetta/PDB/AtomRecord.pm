package Rosetta::PDB::AtomRecord;

use Exporter;

@ISA = qw/ Exporter /;
@EXPORT_OK = qw/
	records_from_file
	records_to_file
	read_first_model
/;


sub records_from_file {
	my $fn = shift;

	my @records;
	open FILE, "<$fn" or die $!;
	while ( my $line = <FILE> ) {
		if ( $line =~ /^ATOM/ ) {
			chomp $line;
			push @records, Rosetta::PDB::AtomRecord->new_from_line($line);
		}
	}

	close FILE or die $!;
	if ( wantarray ) {
		return @records;
	} else {
		return \@records;
	}
}

sub read_first_model {
	my $fn = shift;

	my @records;
	open FILE, "<$fn" or die $!;
	while ( my $line = <FILE> ) {
		if ( $line =~ /^ATOM/ ) {
			chomp $line;
			push @records, Rosetta::PDB::AtomRecord->new_from_line($line);
		} elsif ( $line =~ /ENDMDL/ || $line =~ /TER/ ) {
			last;
		}
	}

	close FILE or die $!;
	if ( wantarray ) {
		return @records;
	} else {
		return \@records;
	}
}

sub records_to_file {
	my $fn      = shift;
	my $records = shift;

	if ( -f $fn ) {
		die "Error: not overwriting file $fn!\n";
	}
	open FILE, ">$fn" or die $!;

	my $last_chain = $records->[0]{chain};
	foreach my $rec (@$records) {
		if ( $last_chain ne $rec->{chain} ) {
			print FILE "TER\n";
		}

		$last_chain = $rec->{chain};
		print FILE $rec->to_string, "\n";
	}
	close FILE or die $!;
}

sub new_from_line {
	my $class = shift;
	my $line  = shift;

	my $self = {};
	bless $self, $class;
	$self->{record_name} = substr( $line,  0, 6 );
	$self->{atom_no}     = substr( $line,  6, 5 );
	$self->{atom_name}   = substr( $line, 12, 4 );
	$self->{res_name}    = substr( $line, 17, 3 );
	$self->{chain}       = substr( $line, 21, 1 );
	$self->{res_no}      = substr( $line, 22, 4 );
	$self->{x}           = substr( $line, 30, 7 );
	$self->{y}           = substr( $line, 38, 7 );
	$self->{z}           = substr( $line, 46, 7 );
	$self->{occupancy}   = substr( $line, 54, 6 );
	$self->{b_factor}    = substr( $line, 60, 6 );

	trim_whitespace_from_entries;

	return $self;
}

sub clone {
	my $self = shift;

	my $copy = {};
	foreach my $key ( keys %$self ) {
		$copy->{$key} = $self->{$key};
	}
	bless $copy;

	return $copy;
}

sub to_string {
	my $self = shift;

	my %format = (
		record_name => '%6s',
		atom_no     => '%5d',
		atom_name   => '%4s',
		res_name    => '%3s',
		res_no      => '%4d',
		chain       => '%1s',
		x           => '%8.3f',
		y           => '%8.3f',
		z           => '%8.3f',
		occupancy   => '%6.2f',
		b_factor    => '%6.2f',
		space       => '%1s',
	);

	my @order = qw/
		record_name atom_no space atom_name space res_name space chain res_no
		space space space space x y z occupancy b_factor
	/;

	my $desired_len = 80;
	my $line = '';
	$self->{space} = '';
	foreach my $key (@order) {
		$line .= sprintf( $format{$key}, $self->{$key} );
	}
	delete $self->{space};

	while ( length($line) < $desired_len ) {
		$line .= ' ';
	}

	return $line;
}

sub trim_whitespace_from_entries {
	my $self = shift;

	for my $k ( keys %$self ) {
		$self->{$k} =~ s/^\s+//;
		$self->{$k} =~ s/\s+$//;
	}

	return $self;
}

sub trim_whitespace {
	my $str = shift;
	$str =~ s/^\s+//g;
	$str =~ s/\s+$//g;
	return $str;
}

1;

__END__

file format definition, taken from:
http://deposit.rcsb.org/adit/docs/pdb_atom_format.html.
 1 -  6        Record name     "ATOM  "
 7 - 11        Integer         Atom serial number.
13 - 16        Atom            Atom name.
17             Character       Alternate location indicator.
18 - 20        Residue name    Residue name.
22             Character       Chain identifier.
23 - 26        Integer         Residue sequence number.
27             AChar           Code for insertion of residues.
31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)       Occupancy.
61 - 66        Real(6.2)       Temperature factor (Default = 0.0).
73 - 76        LString(4)      Segment identifier, left-justified.
77 - 78        LString(2)      Element symbol, right-justified.
79 - 80        LString(2)      Charge on the atom.
