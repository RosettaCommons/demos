package CASP;

sub casp_to_target_id {
	my $id = shift;
	if ( $id =~ /T0(\d{3})/ ) {
		my $target_num = $1;
		return join '', ( 't', $target_num, '_' );
	}
	die "Error: don't recognize id $id!\n";
}

sub parse_domain_id {
	my $id = shift;

	my $parsed_id;
	if ( $id =~ /(t\d+_)(d\d+)/ ) {
		my $target_id = $1;
		my $domain_id = $2;
		$parsed_id = "$target_id$domain_id";
	} elsif ( $id =~ /(t\d+)/ ) {
		my $target_id = $1 . '_';
		my $domain_id = 'd1';

		$parsed_id = "$target_id$domain_id";
	} elsif ( $id =~ /T0(\d+)-(D\d+)/ ) {
		my $target_id = 't' . $1 . '_';
		my $domain_id = lc $2;

		$parsed_id = "$target_id$domain_id";
	} elsif ( $id =~ /T0(\d+)_(\d+)$/ ) {
		my $target_id = 't' . $1 . '_';
		my $domain_id = 'd' . $2;

		$parsed_id = "$target_id$domain_id";
	} elsif ( $id =~ /T0(\d+)$/ ) {
		my $target_id = 't' . $1 . '_';
		my $domain_id = 'd1';
		$parsed_id = "$target_id$domain_id";
	} else {
		#die "Error: don't recognize $id!\n";
		warn "Error: don't recognize $id!\n";
	}

	return $parsed_id;
}

1;
