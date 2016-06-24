FlexPepDock AbInitio Protocol Capture
=====================================

KEYWORDS: DOCKING PEPTIDES

Written by: Barak Raveh, Nir London, Lior Zimmerman, Ora Schueler-Furman  

---

This is a protocol for the de-novo folding and docking of peptides to proteins (a major extension of the previous FlexPepDock refinement protocol)
For this protocol, no initial information of the peptide backbone is requiered, only a placement of an arbitrary peptide strating 
structure of the peptide within the approximate binding pocket.

Setting up the demo
-------------------

Before running the demo, make sure to change the following variables to you 
local environment:

    FlexPepDockAbInitio/prepack_example:                 set PATH_TO_EXE and PATH_TO_DB
    FlexPepDockAbInitio/run_example:                     set PATH_TO_EXE and PATH_TO_DB
    FlexPepDockAbInitio/scripts/frags/make.sh:set        set PATH_TO_EXE and PATH_TO_DB
    FlexPepDockAbInitio/scripts/clustering/cluster.sh:   set PATH_TO_EXE and PATH_TO_DB 
    FlexPepDockAbInitio/scripts/prep_abinitio.sh:        set pathToDemo and pathToVall
    FlexPepDockAbInitio/scripts/frags/make_fragments.pl: set $scratch and $pathToDemo

Running the demo
----------------

These flags are all found in input_files/flags:

IO flags:

    -s start.ppk.pdb                                        # The start structure of the peptide-protein complex
    -native native.pdb                                      # A reference structure for RMSD calculations - MANDATORY!!!
    -out:pdb_gz                                             # silent output flags
    -out:file:silent_struct_type binary
    -out:file:silent decoys.silent
    -scorefile score.sc                                     # name of scorefile

If using multiple processes and no silent file:

    -multiple_processes_writing_to_one_directory

Number of structures to produce (for demo):

    -nstruct 5                                              # number of structures to produce 

Number of structures to produce (for production run):

    -nstruct 50000

FlexPepDock flags:

    -flexPepDocking:lowres_abinitio
    -flexPepDocking:pep_refine                              # Refine after ab-initio
    -flexPepDocking:flexpep_score_only                      # add aditional interesting scores to scorefile

Packing flags:

    -ex1
    -ex2aro
    -use_input_sc
    -unboundrot native.pdb

Fragment picker flags:

    -frag3 frags/frags.3mers.offset
    -frag9 frags/frags.9mers.offset
    -flexPepDocking:frag5 frags/frags.5mers.offset
    -flexPepDocking:frag5_weight 0.25
    -flexPepDocking:frag9_weight 0.1

Example Rosetta Command Line:

    $PATH_TO_EXE/FlexPepDocking.release -database $PATH_TO_DB @flags

Overall protocol execution (demo):

1.  scripts/prep_abinitio.sh 2b1z (preparation step)

    This will create a 'frags' dir with frgments of the peptide for your run. 
    This will also create links to the native.pdb, start.pdb and flags files 
    from the input_files dir.

2.  prepack_example (prepacking step)

    This will create a start.ppk.pdb structure which is the pre-packed complex 
    to start the simulation from. Also, a ppk.score.sc score file for the 
    repacked structure, as well as a prepack.log log of the run.

3.  run_example (docking step)

    This will create the models of the interaction in a silent (compressed) 
    file (5 for the demo, use 50,000 for real life problems) as well as an 
    initial score file for the models. 

4.  scripts/scoring/rescore.sh score.sc (rescoring step)

    This step is NOT needed if you use Rosetta version 3.3 onwards. If you use 
    version 3.2, use this step to create a new score file (named 'newscore.sc') 
    which includes the fpdock-abinitio recommended score for ranking and 
    clustering of the models (reweighted_sc) as the last column. In this case, 
    use 'newscore.sc' instead of 'score.sc' in the next step (clustering).

5.  scripts/clustering/cluster.sh [ntop] 2 score.sc native.pdb decoys.silent reweighted_sc (clustering step)

    This will cluster the [ntop] lowest energy models (use 5 in this demo, 500 
    for real life problems). The script relies on the specified score file 
    'score.sc' (use 'newscore.sc' if using Rosetta version 3.2), and ranks 
    models according to the specified column 'reweighted_sc' (a reweighted 
    version of the Rosetta score), with the specified clustering radius of 2A. 
    The reprasentative models will be shown in the file 
    'clusters_by_reweighted_sc.txt', sorted by the score of their lowest-energy 
    representatives.

Version
-------
Latest version applies to svn revision 45531 (Oct 2011)


References
----------
Raveh B, London N , Zimmerman L & Schueler-Furman O (2011)
Rosetta FlexPepDock ab-initio: Simultaneous Folding, Docking and Refinement of 
Peptides onto their Receptors PLoS ONE 6(4): e18934.
