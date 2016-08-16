#!/usr/bin/perl
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
#  @author James Thompson
#  @author Ray Yu-Ruei Wang, wangyr@u.washington.edu
#

use strict;
use warnings;

use Getopt::Long;

my %options;
$options{renumber} = 1;
$options{outfile}  = 'combined.out';
&GetOptions(
	\%options,
	"outfile=s",
	"renumber!",
);

my @files = @ARGV;

my $tag_generator  = &get_tag_generator;
if ( $options{renumber} ) {
	$tag_generator = &get_renumber_tag_generator;
}

my $printed_header = 0;
open OUT, ">$options{outfile}" or die $!;
foreach my $fn (@files) {
	open FILE, "<$fn" or die $!;
	print $fn;
	my $seen_header = 0;
	my $header = '';
	while ( !$seen_header && (my $line = <FILE>) ) {
		$header .= $line;
		if ( $line =~ /^SCORE:/ ) {
			$seen_header = 1;
			if ( !$printed_header ) {
				$printed_header = 1;
				print OUT $header;
			}
		}
	}

	my $current_tag = 0;
	while ( my $line = <FILE> ) {
		if ( $line =~ /^REMARK/ ) {
			print OUT $line;
		} else {
			chomp $line;
			#my @e   = split /\s+/, $line;
			if ( $line !~ /^(.*)\s+([\/\-\d\w_\.]+)$/ ) {
				warn "Error: don't recognize line $line!\n";
			} else {
				my $lhs = $1;
				my $tag = $2;
				if ( $line =~ /^SCORE/ ) {
					$current_tag = $tag_generator->($tag);
				}
				my $new_line = join ' ', ($lhs,$current_tag);
				print OUT $new_line, "\n";
			}
		}
	}
	close FILE or die $!;
}
close OUT or die $!;

sub get_tag_generator {
	my %tags;
	my $fn = sub {
		my $tag = shift;
		if ( exists $tags{$tag} ) {
			my ($lhs,$num) = split /_/, $tag;

			if ( $tag =~ /([\d\w]+)_([\d+])/ ) {
				$lhs = $1;
				$num = $2;
			}

			if ( !defined $lhs ) {
				print "tag = $tag\n";
				die "Error: lhs $lhs not defined!\n";
			}
			if ( !defined $num ) {
				print "tag = $tag\n";
				die "Error: num $num not defined!\n";
			}

			my $new_tag = join '_', ( $lhs, $num );
			while ( exists $tags{$new_tag} ) {
				$num++;
				$new_tag = join '_', ( $lhs, $num );
			}
			$tags{$new_tag} = 1;
			return $new_tag;
		} else {
			$tags{$tag} = 1;
			return $tag;
		}
	};

	return $fn;
}

sub get_renumber_tag_generator {
	my %tags;
	my $num = 1;
	my $fn = sub {
		my $tag = shift;

		my $lhs = 'S';
		my $new_tag = join '_', ( $lhs, $num );
		while ( exists $tags{$new_tag} ) {
			$num++;
			$new_tag = join '_', ( $lhs, $num );
		}
		$tags{$new_tag} = 1;
		return $new_tag;
	};

	return $fn;
}
