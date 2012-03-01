package Rosetta::Cartesian;

sub cartesian {
	my $first_set = shift @_;
	my @product = map { [ $_ ] } @$first_set;

	foreach my $set (@_) {
		my @new_product;
		foreach my $s (@$set) {
			foreach my $list (@product) {
				push @new_product, [ @$list, $s ];
			}
		}

		@product = @new_product;
	}

	return \@product;
}

1;
