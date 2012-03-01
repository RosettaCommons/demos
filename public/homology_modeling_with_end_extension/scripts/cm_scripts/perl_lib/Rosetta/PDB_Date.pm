package Rosetta::PDB_Date;

#use constant ftp_base => "ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/";
#use constant ftp_link => ftp_base . 'entries.idx';

use File::Basename qw/ dirname /;
my $dir = dirname(__FILE__);
my $dates_file = "$dir/date_map.txt";
#Use date_exclusion_scripts to make date_map.txt file from several entries.idx files 

sub new {
	my $class = shift;

	my $self = {};
	bless $self, $class;
	$self->init( $dates_file );
	return $self;
}

sub all_ids {
	my $self = shift;
	my @ids  = keys %{ $self->{pdb_dates} };

	return \@ids;
}

sub date {
	my $self = shift;
	my $id   = shift;
	$id = lc substr( $id, 0, 4 );

	if ( exists $self->{pdb_dates}{$id} ) {
	    return $self->{pdb_dates}{$id};
	} else {
		return 0;
	}
}

sub init {
	my $self = shift;
	my $file = shift;

	if ( ! -f $file ) {
		$self->update_entries;
	}

	open FILE, "<$file" or die "Error opening file $file ($!)";

	my %pdb_dates;
	while ( my $line = <FILE> ) {
		my ($pdbid,$date) = split /\s+/, $line;
		my ($month, $day, $year);
		if ( $line =~ /\s+(\d+)-(\d+)-(\d+)/ ) {
		    $year  = $1;
		    $month = $2;
		    $day  = $3;
		} else {
		    warn "Error: no date in $line!\n";
		    next;
		}
		$pdbid = lc $pdbid;
		$pdb_dates{$pdbid} = $year . $month . $day;
	}
	close FILE or die $!;
	$self->{pdb_dates} = \%pdb_dates;
}

sub update_entries {
	my $self = shift;
	my $lockfile = $dates_file . '.lock';
	if ( -f $lockfile ) {
		die "Error: not updating: file $lockfile exists!\n";
	}

	use File::Basename qw/ dirname /;

	my $dir = dirname $dates_file;
	unlink $dates_file;
	my $cmd = "cd $dir; wget " . ftp_link;
	print $cmd, "\n";
	system( $cmd );
}

1;
