###
###
### This file is part of the CS-Rosetta Toolbox and is made available under
### GNU General Public License
### Copyright (C) 2011-2012 Oliver Lange
### email: oliver.lange@tum.de
### web: www.csrosetta.org
###
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program.  If not, see <http://www.gnu.org/licenses/>.
###
###
#!/usr/bin/perl

use strict;
use warnings;

# pick_fragments.pl - A program for picking fragments using NMR chemical
# shift data.

use Cwd qw/ abs_path getcwd /;
use File::Path qw/ mkpath /;
use File::Basename;
use Getopt::Long;

use Cwd qw(abs_path);

my %options;
$options{fasta}         = 0;
$options{talos_fn}      = 0;
$options{dir}           = "$ENV{PWD}";
$options{n_procs}       = 4;
$options{frag_sizes} =[];
$options{n_frags}       = 200;
$options{talos_version} = "";
# set these paths now in separate file: setup_paths.pl
# the goal is to keep the pick_fragments.pl clean from any user-specific changes...

#$options{vall}          = "$ENV{HOME}/csrosetta3/frag_picker/vall.dat.2008.apr24.vCS";
#$options{mini_db}       = "$ENV{HOME}/csrosetta3/minirosetta_database";
#$options{blastpgp}      = "$ENV{HOME}/csrosetta3/frag_picker/blast/bin/blastpgp";
#$options{blast_db}      = "$ENV{HOME}/csrosetta3/frag_picker/genomes/nr";
#$options{vall_blast_db} = "$ENV{HOME}/csrosetta3/frag_picker/vall.blast.2008.apr24";
#$options{picker}        = "$ENV{HOME}/csrosetta3/mini/bin/fragment_picker.default.linuxgccrelease";

my $prog_path = dirname(abs_path($0));

require "$prog_path/setup_paths.pl";
&add_paths(\%options);

&GetOptions(
	\%options,
	"fasta=s",
	"talos_fn=s",
	"dir=s",
	"vall=s",
	"mini_db=s",
	"blastpgp=s",
	"frag_sizes=i@",
	"n_procs=i",
	"n_frags=i",
	"talos_version=s",
	"hom!",
);

if ($options{hom}) {
  print "use homologues.\n";
} else {
  print "exclude homologues.\n";
}

#print join ' ', @{$options{frag_sizes}};

logger( "Picking fragments with following options:" );
logger( &options_to_str(\%options) );

# step0 - set up input files in working directory
mkpath($options{dir});
assert_file_exists( "fasta",    $options{fasta}    );
assert_file_exists( "talos_fn", $options{talos_fn} );
assert_dir_exists ( "dir",      $options{dir}      );
copy_files_to_dir( $options{fasta}, $options{talos_fn}, $options{dir} );

my $pred_filename = "$options{fasta}.talos/predSS.tab";
if ( $options{talos_version} eq "pre_2012" ) {
  print "Expecting TALOS VERSION: pre_2012 \n";
  $pred_filename = "$options{fasta}.talos/pred.ss.tab";
} else {
  print "Expecting TALOS VERSION: post_2012 \n";
}

# step1 - run talos
run_talos( $options{talos_fn}, $options{fasta}, $options{dir}, $pred_filename );

# step2 - generate .checkpoint and .homologs file (can be done in parallel with
# step1)
run_blast( $options{fasta}, $options{dir} );

# step3 - pick fragments using minirosetta fragment picker, using inputs from
# step1 and step2.
pick_fragments( $options{fasta}, $options{talos_fn}, $options{dir}, $pred_filename );

##############################
### Subroutine Definitions ###
##############################

sub pick_fragments {
	my $fasta         = shift;
	my $talos_fn      = shift;
	my $dir           = shift;
	my $pred_filename = shift;
	my $picker = $options{picker};

	my $orig_dir = getcwd;
	chdir($dir);

	my $frag_sizes = join ' ', @{$options{frag_sizes}};

	my $flags = <<FLAGS;
#-mute all
-out:level 300
#-out:levels basic.io.database:warning
-mute_info basic.io.database
-database                 $options{mini_db}
-in:file:vall             $options{vall}
-frags:n_frags            $options{n_frags}
-frags:frag_sizes         $frag_sizes
-frags:describe_fragments frags.fsc.score
#-out:file:frag_prefix     $fasta.frags
-out:file:frag_prefix     frags.score
-frags:scoring:config     $fasta.scores.cfg
-in:file:checkpoint       $fasta.checkpoint
-in:file:talos_cs         $talos_fn
-in:file:fasta            $fasta
-frags:ss_pred            $pred_filename talos
-in:file:talos_phi_psi    $fasta.talos/pred.tab
-frags:sigmoid_cs_A       2
-frags:sigmoid_cs_B       4
FLAGS


       open FILE, ">$fasta.pick_flags" or die $!;
	print FILE $flags;
	close FILE or die $!;

	if ( !$options{hom} ) {
	  print "exclude homologues\n";
	  my $flags = <<FLAGS;
-frags:denied_pdb         $fasta.homologs
FLAGS
	  open FILE, ">>$fasta.pick_flags" or die $!;
	  print FILE $flags;
	  close FILE or die $!;

	}

	my $scores_cfg = <<SCORES;
# score name       priority    wght   min_allowed  extras
CSScore             400        1            -
ProfileScoreL1      300        1.5            -
TalosSSSimilarity   200        0.25           -    talos
RamaScore           100        1              -    talos
PhiPsiSquareWell     50        0.15           -
SCORES

	open FILE, ">$fasta.scores.cfg" or die $!;
	print FILE $scores_cfg;
	close FILE or die $!;

	my $cmd = "$picker \@$fasta.pick_flags";
	my $out_fn = join '', ( $fasta, '.frags.', $options{frag_sizes}[0], "mers" );
	produce_output_with_cmd( $cmd, "frags.fsc.score.200.9mers" );

	chdir($orig_dir);
}

sub run_talos {
	my $talos_fn = shift;
	my $fasta = shift;
	my $dir      = shift;
	my $pred_filename = shift;
	my $talos_script = <<TALOS;
mkdir -p $fasta.talos
cd $fasta.talos
talos+ -in ../$talos_fn 2>/dev/null
TALOS

	use Cwd qw/ getcwd /;
	my $orig_dir = getcwd;
	chdir($dir);

	my $talos_script_fn = "talos.script";
	open FILE, ">$talos_script_fn" or die $!;
	print FILE $talos_script;
	close FILE or die $!;

	my $cmd = "bash $talos_script_fn";

	produce_output_with_cmd($cmd,$pred_filename);

	chdir($orig_dir);
}

sub run_blast {
	my $fasta = shift;
	my $dir   = shift;

	use Cwd qw/ getcwd /;
	my $orig_dir = getcwd;
	chdir($dir);

	# run BLAST to make .checkpoint file
	my $chk_file = "$fasta.chk";
	my $cmd = "$options{blastpgp} -t T -i $fasta -F F -j2 -o $fasta.blast -d $options{blast_db} -v10000 -b10000 -K1000 -h0.0009 -e0.0009 -C $chk_file -Q $fasta.pssm -a$options{n_procs}";
	produce_output_with_cmd($cmd,$chk_file);
	my $sequence = read_sequence( $fasta );
	my @matrix   = parse_checkpoint_file( $chk_file );
	@matrix      = finish_checkpoint_matrix( $sequence, \@matrix );
	write_checkpoint_file( "$fasta.checkpoint", $sequence, @matrix );

	# run BLAST to make .homologs file
	my $hom_cmd = "$options{blastpgp} -t T -F F -i $fasta -j 1 -R $chk_file -o $fasta.outn -e 0.05 -d $options{vall_blast_db} -a$options{n_procs}";
	produce_output_with_cmd($hom_cmd,"$fasta.outn");
	my @ids = unique(
		pdb_ids_from_outn("$fasta.outn"),
		pdb_ids_from_blast("$fasta.blast")
	);

	my $query_id = basename($fasta,'.fasta');
	open FILE, ">$fasta.homologs" or die $!;
	foreach my $id (@ids) {
		print FILE join ' ', ( $query_id, $id );
		print FILE "\n";
	}
	close FILE or die $!;

	chdir($orig_dir);
}

sub unique {
	my %h;
	foreach my $d (@_) {
		$h{$d}++;
	}

	my @u = keys %h;
	if ( wantarray ) {
		return @u;
	}
	return \@u;
}


sub pdb_ids_from_blast {
	my $fn = shift;

	my @ids;
	open FILE, "<$fn" or die "Error opening file $fn ($!)";
	while ( my $line = <FILE> ) {
		chomp $fn;
		if ( $line =~ /^>pdb\|([\d\w]{4})\|([\d\w]{1})\s+/ ) {
			my $pdbid = lc $1;
			my $chain = uc $2;
			push @ids, "$pdbid$chain";
		}
	}
	close FILE or die $!;

	if ( wantarray ) {
		return @ids;
	}
	return \@ids;
}

#Vall Blast Database Format
sub pdb_ids_from_outn {
	my $fn = shift;

	my @ids;
	open FILE, "<$fn" or die "Error opening file $fn ($!)";
	while ( my $line = <FILE> ) {
		chomp $fn;
		my @words = split(/ +/, $line);
		next unless (scalar(@words)>= 3);
		if ( $words[2] eq "res" ) {
		    push @ids, $words[0];
		}
	}
	close FILE or die $!;

	if ( wantarray ) {
		return @ids;
	}
	return \@ids;
}

sub produce_output_with_cmd {
	my $cmd = shift;
	my $output_fn = shift;

	logger( "Command: $cmd" );
	if ( -f $output_fn ) {
		logger( "Skipping command, $output_fn exists!" );
		logger( "" );
		return;
	}

#	my $output = `$cmd`;
	system("$cmd");
#	logger( "Output: $output" );
	logger( "Finished running command: $cmd" );

	if ( ! -f $output_fn ) {
		print "Error: expected creation of $output_fn after running '$cmd'!\n";
	}
}

sub copy_files_to_dir {
	my $dir = pop @_;
	use Cwd qw/ abs_path /;
	use File::Copy;

	if ( abs_path($dir) eq abs_path($ENV{PWD}) ) {
		return;
	}

	foreach my $fn (@_) {
		copy( $fn, "$dir/$fn" );

		if ( ! -f "$dir/$fn" ) {
			die "Error: failed to copy $fn to $dir!\n";
		}
	}
}

sub options_to_str {
	my $options = shift;

	my $str = '';

	use List::Util qw/ max /;
	my $max_width  = max map { length($_) } keys %$options;
	my $format_str = "%" . $max_width . "s";
	foreach my $key ( sort keys %$options ) {
		my $val_str = $options->{$key};
		if ( ref($options->{$key}) eq 'ARRAY' ) {
			$val_str = join ' ', @{ $options->{$key} };
		} elsif ( ref($options->{$key}) eq 'HASH' ) {
			my @keys = keys %{ $options->{$key} };
			$val_str = join ', ',
				map { join ': ', ( $_, $options->{$key}{$_} ) } @keys;
		}
		$str .= join ': ', ( sprintf( $format_str, $key ), $val_str );
		$str .= "\n";
	}
	$str .= '-' x 80 . "\n";
	return $str;
}

sub logger {
	for my $str (@_) {
		print $str, "\n";
	}
}

sub assert_file_exists {
	my $name = shift;
	my $fn   = shift;
	if ( ! -f $fn ) {
		die "Error: file doesn't exist! (name = $name, fn = $fn)\n";
	}
}

sub assert_dir_exists {
	my $name = shift;
	my $dir  = shift;
	if ( ! -d $dir ) {
		die "Error: directory doesn't exist! (name = $name, dir = $dir)\n";
	}
}

sub read_sequence {
	my $filename = shift;
	open FILE, "<$filename" or die $!;
	my @file = <FILE>;
	close FILE or die $!;

	my @lines = grep { !/^>/ } map { chomp $_; $_ } @file;
	return join '', @lines;
}

## parse_checkpoint_file -- parses a PSI-BLAST binary checkpoint file.
#
# args:	filename of checkpoint file.
# rets:	N x 20 array containing checkpoint weight values, where N
#				is the size of the protein that BLAST thought it saw...

sub parse_checkpoint_file {
	my $filename = shift;
	my $buf;
	my $seqlen;
	my $seqstr;
	my $i;
	my $j;
	my @aa_order = split(//,'ACDEFGHIKLMNPQRSTVWY');
	my @altschul_mapping = (0,4,3,6,13,7,8,9,11,10,12,2,14,5,1,15,16,19,17,18);
	my @w;
	my @output;

	open(INPUT, $filename) || die ("Couldn't open $filename for reading.\n");

	read(INPUT, $buf, 4) || die ("Couldn't read $filename!\n");
	$seqlen = unpack("i", $buf);

	read(INPUT, $buf, $seqlen) || die ("Premature end: $filename.\n");
	$seqstr = unpack("a$seqlen", $buf);

	for ($i = 0; $i < $seqlen; ++$i) {
		read(INPUT, $buf, 160) || die("Premature end: $filename, line: $i\n");
		@w = unpack("d20", $buf);

		for ($j = 0; $j < 20; ++$j) {
			$output[$i][$j] = $w[$altschul_mapping[$j]];
		}
	}

	return @output;
}

## finish_checkpoint_matrix -- "finishes" the parsed PSI-BLAST checkpoint matrix,
##														 by adding pseudo-counts to any empty columns.
#
# args:	1) sequence string	2) array returned by parse_checkpoint_file
# rets:	"finished" array.	suitable for printing, framing, etc.

sub finish_checkpoint_matrix {
	my $s = shift;
	my $matrix_ref = shift;
	my @matrix = @{$matrix_ref};
	my @sequence = split(//,$s);
	my $i;
	my $j;
	my $line;
	my $sum;
	my @words;
	my @b62;
	my @blos_aa = (0,14,11,2,1,13,3,5,6,7,9,8,10,4,12,15,16,18,19,17);
	my %aaNum = (
		A => 0,
		C => 1,
		D => 2,
		E => 3,
		F => 4,
		G => 5,
		H => 6,
		I => 7,
		K => 8,
		L => 9,
		M => 10,
		N => 11,
		P => 12,
		Q => 13,
		R => 14,
		S => 15,
		T => 16,
		V => 17,
		W => 18,
		Y => 19,
		X => 0	 ### cheap fix for now ##
	);

	(length($s) == scalar(@matrix))
		|| die ("Length mismatch between sequence and checkpoint file!\n");

	my $BLOSUM_PATH = dirname($options{vall_blast_db});
	open(B62,"$BLOSUM_PATH/blosum62.qij")
		|| die "couldnt find blosum62.qij in $BLOSUM_PATH\n";
	$i = 0;

	# read/build the blosum matrix
	while (<B62>) {
		next if ($_ !~ /^\d/);
		chomp;
		@words = split(/\s/);

		for ($j = 0; $j <= $#words; ++$j) {
			$b62[$blos_aa[$i]][$blos_aa[$j]] = $words[$j];
			$b62[$blos_aa[$j]][$blos_aa[$i]] = $words[$j];
		}

		++$i;
	}

	# normalize the blosum matrix so that each row sums to 1.0
	for ($i = 0; $i < 20; ++$i) {
		$sum = 0.0;

		for ($j = 0; $j < 20; ++$j) {
			$sum += $b62[$i][$j];
		}

		for ($j = 0; $j < 20; ++$j) {
			$b62[$i][$j] = ($b62[$i][$j] / $sum);
		}
	}

	# substitute appropriate blosum row for 0 rows
	for ($i = 0; $i <= $#matrix; ++$i) {
		$sum = 0;

		for ($j = 0; $j < 20; ++$j) {
			$sum += $matrix[$i][$j];
		}

		if ($sum == 0) {
			for ($j = 0; $j < 20; ++$j) {
				$matrix[$i][$j] = $b62[$aaNum{$sequence[$i]}][$j];
			}
		}
	}

	return @matrix;
}

sub write_checkpoint_file {
	my ($filename, $sequence, @matrix) = @_;
	my $row;
	my $col;
	my @seq = split(//,$sequence);

	open(OUTPUT, ">$filename");

	die ("Length mismatch between sequence and checkpoint matrix!\n")
		unless (length($sequence) == scalar(@matrix));

	print OUTPUT scalar(@matrix),"\n";

	for ($row = 0; $row <= $#matrix; ++$row) {
		print OUTPUT "$seq[$row] ";
		for ($col = 0; $col < 20; ++$col) {
			printf OUTPUT "%6.4f ", $matrix[$row][$col];
		}
		print OUTPUT "\n";
	}

	#print OUTPUT "END";

	close OUTPUT or die $!;
}
