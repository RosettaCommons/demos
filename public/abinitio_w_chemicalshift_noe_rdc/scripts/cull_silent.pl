#!/usr/bin/perl



use strict;


if ($#ARGV < 0) {
	print STDERR "USAGE: $0 <silentfile> <conditions>\n";
	exit -1;
}

my $silentfile = shift @ARGV;

# parse conditionals
my @sel_fields = ();
my %sel_idx = ();
my @conditionals = ();
my @cutoffs = ();

foreach my $arg ( @ARGV ) {
	if ($arg =~ /(.*)(<=|<|==|>=|>)(.*)/) {
		push @sel_fields, $1;
		push @conditionals, $2;
		push @cutoffs, $3;
	} else {
		print STDERR "Error parsing >> $arg <<\n";
		exit -1;
	}
}


open( INSILENT , $silentfile ) || die "Cannot open $silentfile\n";
my $outfile = $silentfile;
$outfile =~ s/\.[\w]*$/_SELECT\.silent/;

## lines to parse
my $score_field = -1;
my $tag_field = -1;

## info
my $counter = 0;
my @chosen_tags = ();


## phase 1 ... wrie scorefile
while (my $line = <INSILENT>) {
	chomp $line;
	if ($line =~ /^SCORE:/) {
		if ($line =~ /description/) {
			# find the column numbers of score and edens terms
			my @fields = split ' ', $line;
			foreach my $i (0..$#fields) {
				if ($fields[ $i ] eq "descripion^") {
					$tag_field = $i;   # probably always $#fields ...
					print STDERR "Reading tag as column $i\n";
				} else {
					# scan against select fields
					foreach my $sel (@sel_fields) {
						if ($fields[ $i ] =~ /^$sel/) {
							if (defined $sel_idx{ $sel }) {
								print STDERR "Tag '$sel' matches multiple columns!\n";
								exit -1;
							}
							$sel_idx{ $sel } = $i;
							print STDERR "Reading $sel (".$fields[ $i ].") as column $i\n";
						}
					}
				}
			}

			if ( scalar( @sel_fields ) != scalar( keys %sel_idx ) ) {
				print STDERR "Not all fields found!\n";
				exit( -1 );
			}
		} else {
			my @fields = split ' ', $line;
			my $accept = 1;

			foreach my $i (0..$#sel_fields) {
				if ($conditionals[$i] eq '<=' && $fields[$sel_idx{$sel_fields[$i]}]> $cutoffs[$i]) { $accept = 0; }
				if ($conditionals[$i] eq '<'  && $fields[$sel_idx{$sel_fields[$i]}]>=$cutoffs[$i]) { $accept = 0; }
				if ($conditionals[$i] eq '==' && $fields[$sel_idx{$sel_fields[$i]}]!=$cutoffs[$i]) { $accept = 0; }
				if ($conditionals[$i] eq '>'  && $fields[$sel_idx{$sel_fields[$i]}]<=$cutoffs[$i]) { $accept = 0; }
				if ($conditionals[$i] eq '>=' && $fields[$sel_idx{$sel_fields[$i]}]< $cutoffs[$i]) { $accept = 0; }
			}

			if ($accept == 1) {
				#push @chosen_tags, $fields[ $tag_field ];
				push @chosen_tags, $counter;
			}

			$counter++;
		}
	}
}
my $num_read = $counter;
$counter = 0;

## make pruned silent file
`head -2 $silentfile > $outfile`;
open( OUTSILENT , ">>$outfile" );
seek (INSILENT, 0, 0);
my $printstate = 0;
my $counter = 0;
my $tag="XXX";
while (my $line = <INSILENT>) {
	chomp $line;
	if ($line =~ /^SCORE:.*\s(\S+)\s*$/) {
		if ($line !~ /description/) {
			my $istagged = 0;
			foreach my $tag_idx (@chosen_tags) {
				if ( $tag_idx == $counter ) {
					$istagged = 1;
					last;
				}
			}
			$printstate = $istagged;
			$counter++;
			$tag = $1;
		}
	}
	if ($printstate == 1) {
		$line =~ s/ $tag/ S_$counter/;
		print OUTSILENT $line."\n";
	}
}
print STDERR "Selected ".($#chosen_tags+1)." of $num_read structures.\n";
