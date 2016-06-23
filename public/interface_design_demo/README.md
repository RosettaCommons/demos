Interface Design Demo
=====================

KEYWORDS: INTERFACES GENERAL

Author: Kevin Houlihan (khouli at unc dot edu)

This demo contains scripts to redesign residue sequences at the interface of
existing heterodimers. Protocols for running with PackRotamers, MinMover,
and MinMover with a scorefunction produced by optE are included.

Setup
-----

This demo contains scripts with absolute paths that need to be updated. A
script is included to do this for you. (If you don't want to use it a list of
files with absolute paths that need to be set is at the bottom of this
section.)

Run the configure script. Pass your rosetta directory as an argument if it
isn't in your home directory. After setting your rosetta directory, the script
tries to add a weight file.

Example:

    [khouli@killdevil-login2 demo]$ ./configure 
    Use /nas02/home/k/h/khouli/rosetta as Rosetta directory for scripts? [Y/n]y
    /nas02/home/k/h/khouli/rosetta
    Add optE_inf_premin.wts to /nas02/home/k/h/khouli/rosetta/rosetta_database/scoring/weights/? [Y/n]y

Alternatively:

    [khouli@killdevil-login2 demo]$ ./configure other_dir/rosetta
    other_dir/rosetta
    Add optE_inf_premin.wts to other_dir/rosetta/rosetta_database/scoring/weights/? [Y/n]

This fixes absolute filepaths within the interface_design_demo/ directory.
Adding the optE_inf_premin.wts to your weights folder is necessary for the
protocol in scripts/design/minpac_optE_premin/ to run.

The configure script makes a few assumptions. You may want to change these.
What those assumptions are and where they can be found:

- The relative paths downwards from interface_design_demo/ have not been 
  changed e.g. scripts within the demo can access other scripts within the demo 
  at their expected paths

- Your rosetta path contains:

  `rosetta_database/`: referenced by design scripts, e.g. pacrot_s12p_nomin.sh 
  (run by run.sh)

  `rosetta_source/bin/sequence_recovery.linuxgccrelease`: referenced by 
  scripts/analyze/seqRec.sh (run by infpro.sh)

  `rosetta_source/bin/rosetta_scripts.mpi.linuxgccrelease`: referenced by 
  design scripts referenced by scripts/analyze/seqRec.sh (run by infpro.sh)

- Your mpi system is the same as mine and you want to use the same options 
  referenced in design scripts

This can also be done manually. Files with absolute paths:

- design scripts in /scripts/design/ e.g. 
  `scripts/design/pacrot_s12p_nomin/pacrot_s12p_nomin.sh`

- `scripts/analyze/seqRec.sh`

- `inputs/selected_chains/selected_chains.list`

- `inputs/min_nats/nats.list`

Running Design Protocols
------------------------

From within one of the protocol directories in scripts/design/, run the run.sh
script.

What this script is doing:

This script runs the actual design script (the one with the same name as the
directory) and passes it 3 arguments. The first argument names the job after
the directory (and hence creates an identically named sub-directory with output
pdbs), the second passes selected_chains.list as the list of proteins to use,
and the third passes the name of the .xml file containing the Rosetta Script to
be used.

Example:

    [khouli@killdevil-login2 pacrot_s12p_nomin]$ pwd
    /nas02/home/k/h/khouli/interface_design_demo/scripts/design/pacrot_s12p_nomin
    [khouli@killdevil-login2 pacrot_s12p_nomin]$ ./run.sh 
    Group   (-G) : bkuhlman_pi
    Project (-P) : bkuhlman_pi
    Memory Limit (-M) : 4 GB
    Job <287057> is submitted to queue <day>.

Change Rosetta run parameters within the <name-of-protocol>.sh script file.
Change the protocol described by the Rosetta script in the xml file.
Interface residues are determined by the RestrictToInterfaceVector task
operator in the Rosetta script.

Analysis
--------

Once a design run completes, analyze it by running post_run.sh from the same
directory the run.sh script was executed it. Inside the output directory, this
will generate the files top_interface_score, sequence_recovery.txt, and
submatrix.txt as well as some others of much less interest.

What this script is doing:

This script is a convenient way to execute the script
interface_design_demo/scripts/analyze/infpro.sh with the proper working
directory. `../../../analyze/infpro.sh` from within the output directory is
awkward.

In turn, infpro.sh is runs a script to score interfaces of output pdbs,
identify the best trajectories of each protein based on those scores, and
then runs a script to analyze sequence recovery and to aggregate substitutions
in residues between the inputs and designs.

This script assumes your rosetta path contains:

    rosetta_source/bin/sequence_recovery.linuxgccrelease

This can be changed in scripts/analyze/seqRec.sh.

To run infpro.sh, set the directory that conntains the output pdbs as your
working directory and then run interface_design_demo/scripts/analyze/infpro.sh.

Example:

    [khouli@killdevil-login2 pacrot_s12p_nomin]$ pwd
    /nas02/home/k/h/khouli/interface_design_demo/scripts/design/pacrot_s12p_nomin/pacrot_s12p_nomin
    [khouli@killdevil-login2 pacrot_s12p_nomin]$ ls sequencerecovery.txt
    ls: sequencerecovery.txt: No such file or directory
    [khouli@killdevil-login2 pacrot_s12p_nomin]$ ~/interface_design_demo/scripts/analyze/infpro.sh 
    [khouli@killdevil-login2 pacrot_s12p_nomin]$ ls sequencerecovery.txt 
    sequencerecovery.txt

The outputs worth looking at are sequence_recovery.txt, submatrix.txt,
and top_interface_scores.


Example ouputs
--------------

The directory example_outputs/ contains outputs generated by running the
configure script and then running each protocol as described above. Each
protocol produces one output folder as well as a benchmark.txt file. The
benchmark.txt files have been prefixed with their protocol but otherwise all
outputs are exactly as produced within the scripts/design/name_of_protocol/
directories.

Where to change parameters
--------------------------

* Adding/removing input heterodimer structures:

  In the run.sh file for the design, change the variable $input_pdb_list to the 
  filename of a text file containing a list of pdbs to use as inputs. 
  Alternatively, add pdbs to inputs/selected_chains/ and re-run configure. The 
  new pdbs will be added to inputs/selected_chains/selected_chains.list Either 
  way, if you want sequence recovery to function you'll need to add a 
  corresponding structure to the set used for sequence recovery.

* Adding/removing native structures for sequence recovery:

  In the post_run.sh file for the design, change the assignment of 
  $native_pdb_list to the filename of a text file containing a list of pdbs to 
  use as natives. Alternatively, add pdbs to inputs/min_nats/ and re-run 
  configure. The new pdbs will be added to inputs/min_nats/nats.list. If you 
  want sequence recovery to function you'll need to have added a corresponding 
  structure to the set used for input. (They have to align by simple ls, i.e. 
  cap-sensitive non-numerical alphabetical order). 

* Number of trajectories used for sequence recovery:

  In the run.sh file, change the n_traj assignment at the top.

* Options passed to Rosetta and MPI:

  In a protocol directory, edit name_of_protocol.sh

Troubleshooting
---------------

PMGR_COLLECTIVE ERROR following post_run.sh: Fix the path to the 
sequence_recovery binary in `interface_design_demo/scripts/analyze/seqRec.sh`

Epilogue
--------

Happy Rosetta-ing!
