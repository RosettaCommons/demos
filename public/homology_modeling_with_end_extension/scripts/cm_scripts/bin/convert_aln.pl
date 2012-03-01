#!/usr/bin/perl

use strict;
use warnings;

use FindBin;
use Getopt::Long;

use lib "$FindBin::Bin/../perl_lib";
use alignment;
use Rosetta::align_util;

my %options = qw/
	format_in  grishin
	format_out tex
	renumber 1
	max_templates 0
	max_template_pct_id 1
	unique 1
	query_len 0
/;

&GetOptions(
	\%options,
	"format_in=s",
	"format_out=s",
	"renumber!",
	"max_templates=i",
	"max_template_pct_id=f",
	"query_len=i",
	"unique!",
);

my @combined_alns;

my %seen;
my $count = 1;
INFILE: foreach my $infile (@ARGV) {
	if ( ! $infile || ! -f $infile ) {
		die "Error: infile must exist! (given $infile)\n";
	}
	#print STDERR "reading alns from $infile with format $options{format_in} ... ";
	my $alns = Rosetta::align_util::read_alns($infile,$options{format_in},"");
	#print STDERR "done.\n";
	foreach my $aln (@$alns) {
		my $add = 1;
		if ( $options{unique} ) {
			my $aln_str = aln_str($aln);
			$add        = !$seen{$aln_str}++;
		}
		if ( $options{max_template_pct_id} < 1 ) {
			my $pct_id = $aln->get_percent_identity;
			if ( $options{query_len} ) {
				$pct_id *= $aln->ungapped_query_len / $options{query_len};
			}
			$add = ( $add &&
				$pct_id <= $options{max_template_pct_id}
			);
		}

		if ( $add ) {
			push @combined_alns, $aln;
		}
	}
}

#print STDERR "have ", scalar(@combined_alns), " alignments.\n";

if ( exists $combined_alns[0]->{score} ) {
	#@combined_alns = sort { $b->{score} <=> $a->{score} } @combined_alns;
	foreach my $aln (@combined_alns) {
		$aln->{sel_score} = $aln->get_percent_identity * $aln->{score};
	}
	@combined_alns = sort { $b->{sel_score} <=> $a->{sel_score} } @combined_alns;
} else {
	@combined_alns = sort { $a->get_e_value <=> $b->get_e_value } @combined_alns;
}

my $aln_count = 1;
ALN: foreach my $aln (@combined_alns) {
	if ( $options{renumber} ) {
		$aln->template_name(
			join '_', ( substr($aln->template_name,0,5), $aln_count )
		);
	}
	print Rosetta::align_util::convert_aln($aln,$options{format_out}), "\n";
	if ( $options{max_templates} && $aln_count >= $options{max_templates} ) {
		last ALN;
	}
	$aln_count++;
}

sub aln_str {
	my $aln = shift;

	return join "\n", (
		substr( $aln->get_template_name, 0, 5 ),
		$aln->get_query_alignment,
		$aln->get_template_alignment,
	);
}
