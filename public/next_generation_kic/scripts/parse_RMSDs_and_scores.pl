#!/usr/bin/perl -w

###########################################################################
#
# script to parse Rosetta score and local backbone RMSD after 
# remodeling of a segment with the Rosetta loopmodel.executable
#
###########################################################################

use strict;
use sigtrap;
use diagnostics;
use POSIX; # ceil, floor
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

if ($#ARGV < 0) {
    print "error! you need to provide the result directory that contains all the .log[.gz] files\n";
    exit;
}
my ($res_dir) = @ARGV;
my $debug = 1;
my $gnuplot = 0; ## set to 1 to generate a score-vs-RMSD plot (.eps) -- requires having gnuplot installed

my ($l, $logfile, $rmsd, $t_string, $loop_start, $loop_end, $pos, $best_fa, $score_fa, $core_name,
    @p,@fa_acc, @fa_rej, @score_vs_RMSD_fa,, %top_fa);



my @logfiles = `ls $res_dir/*.log*`;
foreach $logfile (@logfiles) {
  chomp $logfile; 
  $core_name = $logfile;
  $core_name =~ s/.gz$//;
  $core_name =~ s/.log$//;
  $best_fa = $score_fa = undef;

  my $fh = read_possibly_gzipped_file($logfile);
  while ($l = <$fh>) {
    chomp $l;
    if ($l =~ /loop_rms (\S+)/) {
      $best_fa = $1;
    } elsif ($l =~ /total_energy (\S+)/) {
      $score_fa = $1;

      if (defined($best_fa)) {
	my $filename_for_match = $core_name . "*_0001.pdb.gz";
	my $filename = `ls $filename_for_match`;
	chomp $filename;
	my $string = join("\t", $score_fa, $best_fa, $filename);
	push @score_vs_RMSD_fa, $string;
	my $m = max(keys %top_fa);
	if (!defined($m) || scalar keys %top_fa < 10 || $score_fa < $m) {
	  $top_fa{$score_fa} = $string; ## we assume scores are unique... 
	  if (scalar keys %top_fa > 10) { ## delete worst case
	    delete $top_fa{$m};
	  }
	}
      }
    }
  }
  close($fh);
}


my $outfile_core = $res_dir;
$outfile_core =~ s/\W+//g;
my $title_core = $outfile_core;
$title_core =~ s/_/-/g;



## Score-vs-RMSD .dat files and (optionally) plots
my ($m, $max_count, @plot_lines, %plot_info);
## append top10
my $datfile2 = $outfile_core . "_fa_score_vs_RMSD.dat";
push @score_vs_RMSD_fa, join("\n", "\n\n## top10", values %top_fa);
$m = min(keys %top_fa);
push @score_vs_RMSD_fa, join("\n", "\n\n## top1", $top_fa{$m}, "");


write_plain_file($datfile2, @score_vs_RMSD_fa);

push @plot_lines, " \"$datfile2\" using 2:1 pt 6 lc rgb \"red\" title \"full-atom\" axes x1y1",
  "\"$datfile2\" index 1:1 using 2:1 pt 13 ps 2 lc rgb \"orange\" title \"fa top10\" axes x1y1", 
  "\"$datfile2\" index 2:2 using 2:1 pt 12 ps 2 lc rgb \"black\" title \"fa top-scoring\" axes x1y1";

my $plotfile = $outfile_core  . "_score_vs_RMSD.plot";
my $eps_file = $outfile_core  . "_score_vs_RMSD.eps";
my $max_rmsd;
my $info = "set title \"" . $title_core . " final models\"
set xlabel \"local backbone RMSD\"
set xtics out
set xrange [0:*]
set ylabel \"Rosetta Energy Units\"
set ytics out
set ytics nomirror
plot " . join(", \\\n ", @plot_lines) . ";\n";
simple_plot($plotfile, $eps_file, $info) unless (!$gnuplot);







sub write_plain_file {
  my ($f, @data) = @_;
  open(DAT, ">$f") || die "Error writing to $f: $!";
  print DAT join ("\n", @data, "");
  close(DAT);
  return;
}


sub simple_plot {
  my ($plotfile, $eps_file, $info) = @_;
  open(PLOT, ">$plotfile") || die "error writing to $plotfile: $!";
 print PLOT "set terminal postscript eps enhanced color \"Helvetica\" 16
set output \"$eps_file\"
$info
";
  close(PLOT);
  print `gnuplot $plotfile`;

  return;
}



## subroutine that opens a (possibly gzipped) file and returns the file handle
## based on http://alumnus.caltech.edu/~svhwan/prodScript/avoidSystemBackticks.html
sub read_possibly_gzipped_file {
  my ($f, $error_msg) = @_;
  if (!defined($error_msg)) {
    $error_msg = " [no error message specified] ";
  }
  if (!-e($f) && -e($f . ".gz")) {
    $f = $f . ".gz";
  } elsif (!-e($f)) {
    print "error -- file $f does not exist -- $error_msg\n";
    return;
  }

  my $INPUTFILE;
  if ($f =~ /\.gz$/) {	# File ends with ".gz"
    my $pid;
    if (not defined($pid = open($INPUTFILE, "-|"))) {
      die "can't fork: $!";
    }
    if ($pid) {
      # parent process - do nothing
    } else {
      # child process
      system( "gzip", "-f", "-c", "-d", $f);
      exit 0;
    }
  } else {
    open($INPUTFILE, $f)
      or die "couldn't open $f: $!";
  }
  return $INPUTFILE;
}

