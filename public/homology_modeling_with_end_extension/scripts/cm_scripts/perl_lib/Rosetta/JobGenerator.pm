package JobGenerator;

use Rosetta::Job;
use Rosetta::Cartesian;
use Text::Template;

sub new {
	my $class = shift;

	my $self = {};
	while ( @_ ) {
		my $key = shift;
		my $val = shift;
		$self->{$key} = $val;
	}

	bless $self, $class;
	return $self;
}

sub generate_jobs {
	my $self = shift;

	my $job_template = $self->{job_template};
	my %vars = %{ $self->vars };

	# split up job_template into a series of jobs
	my @templates;
	my $base = '';
	my $curr_template = '';
	foreach my $line ( split /\n/, $job_template ) {
		if ( $line =~ /^JOB:/ ) {
			if ( length($curr_template) > 0 ) {
				push @templates, $curr_template;
			}
			$curr_template = join "\n", ( $base, "$line\n" );
		} elsif ( length($curr_template) == 0 ) {
			$base .= "$line\n";
		} else {
			if ( $line ne '' ) {
				$curr_template .= "$line\n";
			}
		}
	}
	if ( length($curr_template) > 0 ) {
		push @templates, $curr_template;
	}

	my $total = 1;
	my @sets;
	my @var_order;
	foreach my $var_name ( keys %vars ) {
		push @var_order, $var_name;
		push @sets, $vars{$var_name};
		$total *= scalar @{ $vars{$var_name} };
	}
	print STDERR "generating $total jobs.\n";
	my @product = @{ Rosetta::Cartesian::cartesian( @sets ) };

	my @jobs;
	foreach my $job_template (@templates) {
		my $text_template = Text::Template->new(
			TYPE => 'STRING', SOURCE => $job_template,
		);

		foreach my $combination( @product ) {
			my %data;
			for my $idx ( 0 .. scalar(@var_order) - 1 ) {
				my $name     = $var_order[$idx];
				my $value    = $combination->[$idx];
				$data{$name} = $value;
			}

			my $job_str = $text_template->fill_in( HASH => \%data );
			my $job = Rosetta::Job->new;
			$job->from_string( $job_str );
			push @jobs, $job;
		}
	}

	return \@jobs;
}

sub job_template {
	my $self = shift;
	if ( @_ ) {
		my $templ = shift;
		$self->{job_template} = $templ;
	}

	return $self->{job_template};
}

sub vars {
	my $self = shift;
	if ( @_ ) {
		my $vars = shift;
		$self->{vars} = $vars;
	}
	return $self->{vars};
}

sub add_var {
	my $self = shift;
	my $var  = shift;
	my $var_list = shift;

	my @copy = @$var_list;

	$self->{vars}{$var} = \@copy;
}

1;
