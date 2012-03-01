package Potential;

sub init_bounds {
	my $self = shift;
	my $bounds_fn = shift;
	open BOUNDS, "<$bounds_fn" or die $!;
	my @file = <BOUNDS>;
	close BOUNDS or die $!;

	my %bins;
	foreach my $line (@file) {
		chomp $line;
		my ($tag,@bounds) = split /\s+/, $line;
		$bins{$tag} = \@bounds;
	}
	$self->{bins} = \%bins;
}

sub bin_order {
	my $self = shift;
	my @bins = @{ $self->bins() };
	my %bin_orders;
	my $idx = 1;
	foreach my $bin (@bins) {
		$bin_orders{$bin} = $idx;
		$idx++;
	}

	return \%bin_orders;
}

sub blur_code {
	my $code = shift;
	my @d = split /\./, $code;
	pop @d;

	return join '.', @d;
}

sub mean {
	my $vals = shift;

	my $total = 0;
	foreach my $d (@$vals) {
		$total += $d;
	}

	my $mean = $total / scalar(@$vals);
	return $mean;
}

sub get_bin_index {
	my $self = shift;
	my $var  = shift;
	my $val  = shift;

	if ( !defined $self->{bins}{$var} ) {
		die "Error: don't know anything about var $var!";
	}

	my $count = 0;
	foreach my $b ( @{ $self->{bins}{$var} } ) {
		if ( $val <= $b ) {
			return $count;
		}
		$count++;
	}

	return $count;
}

1;
