package blast;

use alignment;

my $id_regex = "[\d\w_\|]";

sub parse_alignments {
	my $file = shift;

	if ( -z $file ) {
		die "Error: empty file $file!\n";
	}

	open FILE, "<$file" or die $!;
	my @file = <FILE>;
	close FILE or die $!;

	my @alignments;
	my $aln;
	my $query_name;
	my $total_query_len;
	foreach my $line (@file) {
		chomp $line;
		if ( $line =~ /^\s*$/ ) { next; }
		if ( $line =~ /^Query=\s*($id_regex+)/ ) {
			$query_name = $1;
		}

		if ( $line =~ /^>([\w\d_|]+)\s*/ ) {
			if ( exists $aln->{e_value} ) { push @alignments, $aln; }

			# make a new $aln
			my %nested_aln;
			$aln = \%nested_aln;
			$aln->{template_name}   = $1;
			$aln->{query_name}      = $query_name;
			$aln->{total_query_len} = $total_query_len;

# Score =	362 bits (930), Expect = e-101
		} elsif (
			$line =~ /Score\s*=\s*([\d\.]+) bits \((\d+)\),\s*Expect\s*=\s*([\.e\d-]+)/ ) {
			#$line =~ /Score\s*=\s*([\d\.]+) bits \((\d+)\),\s*Expect\s*=/ )

			$aln->{score}	= $1;
			$aln->{bits}	= $2;
			$aln->{e_value} = $3;
			if ( substr( $aln->{e_value}, 0, 1 ) eq 'e' ) {
				$aln->{e_value} = '1'.$aln->{e_value};
			}
#Probab=100.00	E-value=0	Score=384.97	Aligned_columns=112	Identities=100%
		} elsif ( $line =~ /Query:*\s*(\d+)\s+([\w-]+)\s+(\d+)/ ) {
			my $qt	 = 'Q';
			my $start = $1;
			my $data	 = $2;
			my $stop	 = $3;

			$aln->{$qt}->{alignment} .= $data;
			# define start and stop appropriately
			if ( ! exists $aln->{$qt}->{start} ) {
				$aln->{$qt}->{start} = $start;
			}
			$aln->{$qt}->{stop} = $stop;
			$aln->{$qt}->{len} += length($data);
		} elsif ( $line =~ /Sbjct:*\s*(\d+)\s+([\w-]+)\s+(\d+)/ ) {
			my $qt	 = 'T';
			my $start = $1;
			my $data	 = $2;
			my $stop	 = $3;

			$aln->{$qt}->{alignment} .= $data;
			# define start and stop appropriately
			if ( ! exists $aln->{$qt}->{start} ) {
				$aln->{$qt}->{start} = $start;
			}
			$aln->{$qt}->{stop} = $stop;
			$aln->{$qt}->{len} += length($data);

		#} elsif ( $line =~ /\((\d+)\)\s+letters/ ) {
		} elsif ( $line =~ /\((\d+) letters\)/ ) {
			$total_query_len = $1;
		} else {
			#print "Error: don't recognize line:\n $line\n";
		}
	}

	if ( defined $aln ) {
		push @alignments, $aln;
	}

	my @objects = map { alignment->new( $_ ) } @alignments;
	return \@objects;
}

1;
