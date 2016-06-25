# optE: Refit Rosetta reference energies

KEYWORDS: UTILITIES GENERAL

This demo contains the input files and command lines
necessary to refit the reference energies for the
Talaris2013 score function using optE's sequence-
profile-recovery procedure.

The demo should be run from within the
  run/
directory, using the launch_optE.scr script.  The
demo is set up to run on 27 processors, but can be
run on fewer processors.  If the demo is running
too slowly, comment out the -ex1 and -ex2 flags
from the optE_seqprof.flags file. There  are several
things about the launch_optE.scr script that will have
to be edited.

1) the path to the executable, in launch_optE.scr
2) the mpi launch command, in launch_optE.scr, and
3) the path to the rosetta database, in score.flags

After optE has finished running, use the python script,
find_best_weight_set.py, in the
  scripts/
directory, to determine which of the iterations was
optimal.

$> $ROSETTA3/bin/optE_parallel.default.linuxgccrelease @optE_seqprof.flags.short
