package DistancePotential;

use Potential;
use FindBin qw($Bin);

@ISA = qw/ Potential /;

use constant bounds_file     => $Bin.'/../cm_data/bounds.txt';
use constant table_data_file => $Bin.'/../cm_data/hhsearch_table.txt';

sub new {
	my $class = shift;
	my $self  = {};

	bless $self, $class;

	$self->init_bounds( bounds_file );
	$self->init_table_data( table_data_file );

	return $self;
}

sub gen_pairs {
	my $length = shift;

	use constant seqsep => 10;
	my @pairs;
	for my $resi ( 1 .. $length ) {
	for my $resj ( $resi + 1 .. $length ) {
		if ( abs( $resj - $resi ) >= seqsep ) {
			push @pairs, [ $resi, $resj ];
		}
	}
	}

	return \@pairs;
}

sub init_table_data {
	my $self = shift;
	my $file = shift;

	open TABLE, "<$file" or die $!;
	my @file = <TABLE>;
	close TABLE or die $!;

	my %table_data;
	my $header = shift @file;
	foreach my $line (@file) {
		chomp $line;
		my @d = split /\s+/, $line;
		my $sdev = pop @d;
		pop @d; # remove pred_dist
		my $code = join '.', @d;
		if ( exists $table_data{$code} ) {
			use List::Util qw/ min /;
			$sdev = min( $table_data{$code}, $sdev );
		}
		$table_data{$code} = $sdev;
	}

	$self->{table_data} = \%table_data;
}

sub get_sdev {
	my $self = shift;
	my $v = shift;

	my $code = $self->get_code( $v );
	if ( !exists $self->{table_data}{$code} ) {
		$self->get_best_match( $code, \%table_data );
		#use Data::Dumper::Simple;
		#print Dumper( $v );
		#die "Error no code for $code!\n";
	}

	my $sdev = $self->{table_data}{$code};

	if ( !defined $sdev ) {
		#print "code = $code\n";
		#use Data::Dumper::Simple;
		#print Dumper( $self );
		#exit 1;
		$sdev = 4.0;
	}

	return $sdev;
}

sub get_best_match {
	my $self = shift;
	my $code = shift;
	my $data = shift;

	my $new_code = $self->blur_code( $code );
	my @candidates = grep { $_ =~ /^$new_code/ } keys %$data;
	my @sdevs = map { $data->{$_} } @candidates;

	my $average_sdev = 4.0;
	if ( scalar(@sdevs) != 0 ) {
		$average_sdev = $self->mean( \@sdevs );
	}

	return $average_sdev;
}

sub get_code {
	my $self = shift;
	my $hashref = shift;

	my @codes;
	foreach my $key ( sort keys %$hashref ) {
		my $val = $hashref->{$key};
		my $bin = $self->get_bin_index( $key, $val );
		push @codes, $bin;
	}
	my $code = join ".", @codes;
	return $code;
}

sub weight_predictions {
	my $preds = shift;
}

sub transform_weights {
	my $weights = shift;

	my $total = 0;
	foreach my $w (@$weights) {
		$total += $w;
	}

	my @transformed = map { $_ / $total } @$weights;
	return \@transformed;
}


1;
