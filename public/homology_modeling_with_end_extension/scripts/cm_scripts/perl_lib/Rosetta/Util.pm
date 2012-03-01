package Rosetta::Util;
require Exporter;

use Cwd;
use FindBin;
use File::Copy;
use File::Path;
use File::Basename;

@ISA = qw/ Exporter /;
@EXPORT = qw/
	glob1
	unique
	mkdir_safe
	mkdir_recursive
	copy_safe
	symlink_safe
	assemble_path
	gzip_files
	glob_safe
	basename_safe
	cwd
	make_local_link
	run_cmd
	print_to_file
	read_vars
	file_exists
	script_dir
	cm_base
	eval_template_str
	trim_whitespace
	options_to_str
	load_options_from_file
	assert_files_exist
	find_files_by_regex
/;

sub run_cmd {
	my $cmd = shift;
	my $output = `$cmd 2>&1`;
	return $output;
}

sub make_local_link {
	my $dir = shift;
	my $file1 = shift;
	my $file2 = shift;

	my $orig_dir = getcwd;
	chdir( $dir );

	if ( ! -l $file2 ) {
		symlink_safe( $file1, $file2 );
	}

	if ( ! -l $file2 ) {
		die "Error making link from $file1 -> $file2 in $dir!\n";
	}
	chdir $orig_dir;
	return;
}

sub cwd {
	return getcwd;
}

sub print_to_file {
	my $fn  = shift;
	my $str = shift;

	open FILE, ">$fn" or die "Error opening file $fn ($!)!\n";
	print FILE $str;
	close FILE or die "Error closing file $fn ($!)\n";
}

sub basename_safe {
	my $fn = shift;
	return basename $fn;
}

sub unique {
	my $array = shift;

	my %hash;
	foreach my $a (@$array) {
		$hash{$a}++;
	}

	my @unique = keys %hash;
	if ( wantarray ) {
		return @unique;
	}
	return \@unique;
}

sub mkdir_recursive {
	my $full_dir  = shift;
	mkpath( $full_dir );

	if ( ! -d $full_dir ) {
		die "Error: creating directory $full_dir!\n";
	}
}

sub mkdir_safe {
	my $dirname = shift;

	use File::Path;
	mkpath $dirname;
	if ( ! -d $dirname ) {
		die "Error: can't mkdir $dirname!\n";
	}
}

sub copy_safe {
	my $filename = shift;
	my $dest_dir = shift;

	if ( -e assemble_path( $dest_dir, $filename ) ) {
		#print STDERR "not overwriting ",
		#	assemble_path( $dest_dir, $filename ), "\n";
		return;
	}

	my $error = 0;
	my $errstr = '';
	if ( ! $filename || ! -f $filename ) {
		$error++;
		$errstr .= "Error with file (given $filename)!\n";
	}
	if ( ! $dest_dir || ! -d $dest_dir ) {
		$error++;
		$errstr .= "Error with dir (given $dest_dir)!\n";
	}

	$errstr .= "Error copying filename $filename to $dest_dir!\n";

	if ( $error ) {
		die $errstr;
	}

	use Cwd qw/ abs_path /;
	if ( abs_path($filename) ne abs_path(assemble_path($dest_dir)) ) {
		copy( $filename, $dest_dir ) or die $errstr;
	}
}

sub gzip_files {
	foreach my $fn (@_) {
		# gzip without removing original
		my $gz_file = join '.', ($fn,'gz');
		if ( ! -f $gz_file ) {
			print "gzipping $fn -> $gz_file\n";
			my $cmd = "gzip <$fn > $gz_file";
			system( $cmd );
		}
	}
}

sub assemble_path {
	my $str = join '/', @_;
	return $str;
}

sub symlink_safe {
	my $orig = shift;
	my $link = shift;

	symlink( $orig, $link );
	if ( ! -l $link ) {
		die "Error creating link $link!\n";
	}
}

sub glob_safe {
	my $dir = shift;
	my $pat = shift;

	if ( !$dir || ! -d $dir ) {
		die "Error in glob_safe: no dir named $dir!\n";
	}

	if ( !$pat ) {
		die "Error in glob_safe: no pattern provided!\n";
	}

	my @files = glob( assemble_path( $dir, $pat ) );
	if ( scalar(@files) == 0 ) {
		print "Error: no files found in $dir (searched with $pat)!\n";
	}

	if ( wantarray ) {
		return @files;
	}
	return \@files;
}

sub read_vars {
	my $fn = shift;

	my %dat;
	open FILE, "<$fn" or die $!;
	while ( my $line = <FILE> ) {
		chomp $line;
		my @d = split /\s+/, $line;
		my $name = shift @d;

		$dat{$name} = \@d;
	}

	return \%dat;
}

sub file_exists {
	my $fn = shift;
	return ( defined $fn && -f $fn );
}

sub script_dir {
	use FindBin;
	my $dir = "$FindBin::Bin/bin";
	if ( $dir =~ m/bin\/bin$/ ) {
		$dir = "$FindBin::Bin";
	}

	return $dir;
}

sub cm_base {
	use FindBin;
	my $dir = "$FindBin::Bin";
	if ( basename($dir) eq 'tests' || basename($dir) eq 'bin' ) {
		$dir = dirname($dir);
	}
	return $dir;
}

sub eval_template_str {
	my $str  = shift;
	my $data = shift;

	use Text::Template;
	my $templ = Text::Template->new( TYPE => 'STRING', SOURCE => $str );
	return $templ->fill_in( HASH => $data );
}

sub glob1 {
	my $pat = shift;

	my @fns = glob( $pat );
	return $fns[0];
}

sub trim_whitespace {
	my $str  = shift;
	my $copy = $str;
	$copy =~ s/^\s+//g;
	$copy =~ s/\s+$//g;

	return $copy;
}

sub options_to_str {
	my $options = shift;

	my $str = '';
	foreach my $key ( sort keys %$options ) {
		$str .= "\t";
		my $val_str = $options->{$key};
		if ( ref($options->{$key}) eq 'ARRAY' ) {
			$val_str = join ' ', @{ $options->{$key} };
		}
		$str .= join ': ', ( $key, $val_str );
		$str .= "\n";
	}
	$str .= '-' x 80, "\n";
	return $str;
}

sub load_options_from_file {
	my $fn   = shift;
	my $opts = shift;

	if ( ! -f $fn ) {
		die "Not reading options from $fn as it doesn't exist!\n";
	}

	open FILE, "<$fn" or die $!;
	while ( my $line = <FILE> ) {
		chomp $line;
		if ( $line !~ /^\s*$/ && $line !~ /^\s*#/ ) {
			my ($key,@vals) = split /\s+/, $line;
			@vals = grep { !/^#/ } @vals;

			# check for data type
			if ( exists $opts->{$key} && ref($opts->{$key}) eq 'ARRAY' ) {
				$opts->{$key} = \@vals;
			} else {
				my $val = join ' ', @vals;
				$opts->{$key} = $val;
			}
		}
	}
	close FILE or die $!;

	return $opts;
}

sub assert_files_exist {
	foreach my $fn (@_) {
		if ( ! -f $fn ) {
			die "Error: file $fn doesn't exist!\n";
		}
	}
}

sub find_files_by_regex {
   my $path  = shift;
   my $regex = shift;

	use Cwd qw/ abs_path /;

   my @files;
   my $wanted = sub {
      if ( $_ =~ /$regex/ ) {
         push @files, abs_path($_);
      }
   };
   find( $wanted, $path );

   if ( wantarray ) {
      return @files;
   }
   return \@files;
}

1;
