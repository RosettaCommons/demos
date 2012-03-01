package Rosetta::Protocol;

use Cwd qw/ abs_path /;
use File::Basename;
use Text::Template;

use Rosetta::Job;
use Rosetta::Util;

sub new {
	my $class = shift;
	my $flag_template_fn = shift;

	if ( !defined $flag_template_fn || ! -f $flag_template_fn ) {
		my $errstr = "Error: must provide template file (given $flag_template_fn!)\n";
		if ( $flag_template_fn =~ /\.template\.template/ ) {
			my $temp_name = $flag_template_fn;
			$temp_name =~ s/\.template\.template/\.template/g;
			$errstr .= "(did you mean to provide $templ_name?)\n";
		}
		die $errstr;
	}

	my $self = bless {}, $class;
	$self->flag_template_fn( $flag_template_fn );
	return $self;
}

sub make_flags {
	my $self = shift;
	my $target_files = shift;
	my $template = Text::Template->new(
		TYPE => 'STRING', SOURCE => $self->flag_template()
	);
	my $flag_str = $template->fill_in( HASH => $target_files );
	return $flag_str;
}

# Assumes that whitespace isn't significant, and that File::Basename::fileparse
# can tell which files are which.
sub make_boinc_flags {
	my $self         = shift;
	my $target_files = shift;
	my $flags        = $self->make_flags($target_files);

	# remove anything that looks like a path. Alternatively, we could just
	# parse the paths from $target_files, but people might add files to the
	# templates.
	my @lines = split /\n/, $flags;
	my @boinc_lines;
	foreach my $line (@lines) {
		chomp $line;
		my @tokens = split /\s+/, $line;

		my @new_tokens;
		foreach my $token (@tokens) {
			my ($base,$path,$type) = fileparse( $token );
			push @new_tokens, $base;
		}

		my $boinc_line = join ' ', @new_tokens;
		push @boinc_lines, $boinc_line . "\n";
	}

	my $str = join '', @boinc_lines;
	return $str;
}

sub boinc_tag {
	my $self = shift;
	if ( @_ ) {
		$self->{boinc_tag} = shift;
	}

	return $self->{boinc_tag};
}

sub save_all_out {
	my $self = shift;
	if ( @_ ) {
		$self->{save_all_out} = shift;
	}

	return $self->{save_all_out};
}

sub setup_boinc_job {
	my $self   = shift;
	my $dir    = shift;
	my $files  = shift;

	my $username  = $ENV{USER};
	my $target_id = $files->value( 'target_id' );

	use File::Path;
	mkpath $dir;

	my $method_name = basename $self->flag_template_fn;
	$method_name =~ s/\.template//g;

	my $unique_boinc_name = join '.', (
		join '_', ($target_id,'boinc'),
		$self->boinc_tag,
		$method_name,
		$username,
	);

	my $results_tag = 'IGNORE_THE_REST';
	if ( $self->{save_all_out} ) {
		$results_tag = 'SAVE_ALL_OUT';
	}
	my $unique_boinc_job_name = join '_', ($unique_boinc_name,$results_tag);

	#print STDERR "setting up $unique_boinc_name in $dir\n";

	my $flags_fn = join '.', ($unique_boinc_name,'boinc','flags');
	$flags_fn = assemble_path( $dir, $flags_fn );

	if ( ! -f $flags_fn ) {
		my $str = $self->make_boinc_flags($files);
		open FILE, ">$flags_fn" or die "Error opening file $flags_fn ($!)";
		print FILE $str, "\n";
		close FILE or die "Error closing file $flags_fn ($!)";
	}

	my $zip_fn = join '.', ($unique_boinc_name,'boinc','zip');
	$zip_fn = assemble_path( $dir, $zip_fn );
	if ( ! -f $zip_fn ) {
		$self->make_zip_file( $files, $zip_fn );
	}

	# use relative paths inside of the submit script
	$zip_fn        = abs_path( $zip_fn );
	$flags_fn      = abs_path( $flags_fn );
	$zip_fn_base   = basename( $zip_fn );
	$flags_fn_base = basename( $flags_fn );

	#my $input_file_str = "inputfiles = $flags_fn,$zip_fn";
	my $input_file_str = "inputfiles = " .
		join (',', map { abs_path($_) } ( $zip_fn, $flags_fn ) );

	my $boinc_job_str = <<BOINC_JOB;
application = minirosetta

name = $unique_boinc_job_name
description = Comparative Modeling Benchmark Run: $target_id using protocol $method_name on target $target_id by user $username
$input_file_str
arguments =  \@$flags_fn_base -silent_gz -mute all -out:file:silent default.out -in:file:boinc_wu_zip $zip_fn_base
resultfiles = default.out.gz
queue = 10

BOINC_JOB

	my $boinc_job_fn = assemble_path(
		$dir, join '.', ($unique_boinc_name,'boinc.job')
	);
	if ( ! -f $boinc_job_fn ) {
		open FILE, ">$boinc_job_fn" or die "Error opening file $boinc_job_fn ($!)";
		print FILE $boinc_job_str;
		close FILE or die $!;
	}

	my $files_created = {
		flags_fn => $flags_fn,
		zip_fn   => $zip_fn,
		job_fn   => $boinc_job_fn,
	};

	foreach my $name ( keys %$files_created ) {
		$files_created->{$name} = abs_path( $files_created->{$name} );
	}

	#print STDERR "finished BOINC job setup.\n";

	return $files_created;
}

sub flag_template_fn {
	my $self = shift;
	if ( @_ ) {
		$self->{flag_template_fn} = shift;
	}

	return $self->{flag_template_fn};
}

sub make_zip_file {
	my $self = shift;
	my $target_files = shift;
	my $zipfile_name = shift;

	my $zipfile_dir;
	if ( !defined $zipfile_dir ) {
		$zipfile_dir = $ENV{PWD};
	}

	my $flags = $self->make_flags($target_files);
	my @lines = split /\n/, $flags;

	my @files;
	foreach my $line (@lines) {
		chomp $line;
		my @tokens = split /\s+/, $line;

		foreach my $token (@tokens) {
			if ( -f $token ) {
				push @files, $token;
			}
		}
	}
	if ( exists $self->{extra_files} ) {
		foreach my $fn (@{$self->{extra_files}}) {
			if ( ! (-f $fn || -d $fn) ) {
				die "Error: can't find fn $fn!\n";
			}
			push @files, $fn;
		}
	}

	#print STDERR "making $zipfile_name in $zipfile_dir\n";
	#print STDERR "zipping files into $zipfile_name\n";
	#print STDERR join "\n", map { "\tzip($_)" } @files;
	#print STDERR "\n";

	my $zip_job = Rosetta::Job->new(
		executable => 'zip',
		args => [ '-j', $zipfile_name, @files ],
		lockfile => "$zipfile_name.lock",
		results_file => $zipfile_name,
		dir => $ENV{PWD},
	);
	$zip_job->run;
}

sub flag_template {
	my $self = shift;

	if ( !defined $self->{flag_template} ) {
		my $flag_template_fn = $self->{flag_template_fn};
		open FILE, "<$flag_template_fn" or die $!;
		my @file = <FILE>;
		close FILE or die $!;

		my $str = join '', @file;
		$self->{flag_template} = $str;
	}

	return $self->{flag_template};
}

1;
