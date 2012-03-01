package Statistics;

use Exporter;

@ISA = qw/ Exporter /;
@EXPORT_OK = qw/
	mean
	calc_percentile
/;

sub mean {
	use List::Util qw/ sum /;
	my $sum = sum @_;
	return ( $sum / scalar(@_) );
}

sub calc_percentile {
	my $data       = shift;
	my $real_col   = shift;
	my $label_col  = shift;
	my $percentile = shift;

	if ( !defined $percentile || $percentile < 0 || $percentile > 1 ) {
		die "Error: need a valid percentile (given $percentile)!\n";
	}

	my %pct_dat;
	foreach my $d (@$data) {
		push @{ $pct_dat{$d->{$label_col}} }, $d->{$real_col};
	}

	my $pct_name = join '_', ( 'pct', $percentile, $real_col );
	my @new_data;
	foreach my $key ( keys %pct_dat ) {
		my $N   = scalar( @{ $pct_dat{$key} } );
		my $idx = int( $percentile * $N + 0.5 );
		my @d   = sort { $a <=> $b } @{ $pct_dat{$key} };
		#$pct_dat{$key} = $d[$idx];
		push @new_data, {
			$label_col => $key,
			#real_col   => $real_col,
			N          => $N,
			$pct_name  => $d[$idx],
		}
	}

	return \@new_data;
}

1;
