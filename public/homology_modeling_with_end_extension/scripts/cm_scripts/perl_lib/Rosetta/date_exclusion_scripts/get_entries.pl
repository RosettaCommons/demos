#!/usr/bin/perl

use strict;
use warnings;

use LWP::Simple;

my @snapshot_dates =  qw/
	20050106 20060103 20070102 20070731 20080107 20090105 20090316 20091201 20100104 20110103 20100101
/;

my $base_url = "ftp://snapshots.rcsb.org";

#ftp://snapshots.rcsb.org/20070102/pub/pdb/derived_data/index
foreach my $date (@snapshot_dates) {
	my $url = join '/', ( $base_url, $date, 'pub/pdb/derived_data/index/entries.idx' );
	my $output_fn = "$date.entries.idx";
	if ( ! -f $output_fn ) {
		my $cmd = "wget $url -O $output_fn --quiet";
		my $output = `$cmd`;
	}
}

# download current entries.idx as well
my $output_fn = 'entries.idx';
if ( ! -f $output_fn ) {
	my $url = 'ftp://ftp.wwpdb.org:/pub/pdb/derived_data/index/entries.idx';
	my $cmd = "wget $url -O $output_fn --quiet";
	system($cmd);
}

my %pdb_dates;
my %pdb_date_lines;
my @files = map { "$_.entries.idx" } @snapshot_dates;
push @files, 'entries.idx';
foreach my $fn (@files) {
	open FILE, "<$fn" or die $!;
	<FILE>; <FILE>;
	while ( my $line = <FILE> ) {
		chomp $line;
		my ($pdbid,$header,@rest) = split /\s+/, $line;
		my ($month,$day,$year);
		if ( $line =~ m#\s*(\d{2})/(\d{2})/(\d{2})\s+# ) {
			$month = $1;
			$day   = $2;
         $year  = $3;
      } else {
			print "line = $line\n";
         die "Error: no date in $line!\n";
         next;
      }

      # stupid years!
      if ( $year - 50 < 0 ) {
         $year += 2000;
      } else {
         $year += 1900;
      }
		$pdbid = lc $pdbid;
		if ( exists $pdb_dates{$pdbid} ) {
			my $new_date = join '', ( $year, $month, $day );
			if ( $pdb_dates{$pdbid} != $new_date ) {
				#print "updating date for $pdbid\n";
				#print "\t", $pdb_dates{$pdbid}, " => ", $new_date, "\n";
				#use List::Util qw/ min /;
				$pdb_dates{$pdbid} = min( $new_date, $pdb_dates{$pdbid} );
				if ( $pdb_dates{$pdbid} > $new_date ) {
					$pdb_date_lines{$pdbid} = $line;
				}
			}
		} else {
			$pdb_date_lines{$pdbid} = $line;
		}

		#$pdb_dates{$pdbid} = $year . $month . $day;
	} # line = <FILE>
	close FILE or die $!;
} # @files

open FILE, "<pdb_seqres.txt" or die $!;
my %missing;
while ( my $line = <FILE> ) {
	if ( $line =~ /^>([\d\w]{4})_/ ) {
		chomp $line;
		my $pdbid = lc $1;
		if ( !exists $pdb_dates{$pdbid} ) {
			#print "missing date for $pdbid\n";
			$missing{$pdbid}++;
		}
	} elsif ( $line =~ /^>/ ) {
		die "didn't match $line!\n";
	}
}

foreach my $id ( keys %missing ) {
	print "missing date for $id.\n";
}

open FILE, ">pdb_dates.txt" or die $!;
print FILE join ' ', qw/ pdbid pdb_date /;
print FILE "\n";
foreach my $pdbid ( sort { $pdb_dates{$a} <=> $pdb_dates{$b} } keys %pdb_dates ) {
	print FILE $pdb_date_lines{$pdbid}, "\n";
	#print FILE join ' ', ( $pdbid, $pdb_dates{$pdbid} );
	#print FILE "\n";
}
close FILE or die $!;
