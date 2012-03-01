package Rosetta::AA;

my @amino_acids = qw/
	A C D E F G H I K L M N P Q R S T V W Y
/;

my %map_three_to_one = qw/
	ALA A ARG R ASN N ASP D CYS C GLU E GLN Q GLY G HIS H ILE I
	LEU L LYS K MET M PHE F PRO P SER S THR T TRP W TYR Y VAL V
/;
my %map_one_to_three = reverse %map_three_to_one;

sub name3_to_name1 {
	my $name3 = shift;

	my $name1 = 'X';
	if ( exists $map_three_to_one{$name3} ) {
		$name1 = $map_three_to_one{$name3};
	}

	return $name1;
}

sub is_valid_name1 {
	my $aa = shift;
	if ( exists $map_one_to_three{$aa} ) {
		return 1;
	}
	return 0;
}

sub name1_to_name3 {
	my $name1 = shift;

	my $name3 = 'UNK';
	if ( exists $map_one_to_three{$name1} ) {
		$name1 = $map_one_to_three{$name1};
	}

	return $name3;
}

1;
