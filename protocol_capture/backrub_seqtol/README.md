Backrub Sequence Tolerance
==========================

KEYWORDS: DESIGN GENERAL

Author: Colin A. Smith  
Protocol backrub_seqtol

---

This protocol is designed to predict the tolerated sequence space for a
given protein-protein interface or protein domain. It involves generating an
ensemble of backbone structures by running backrub simulations on an input
structure. For each member of the ensemble, a large number of sequences are 
scored and then Boltzmann weighted to generate a position weight matrix for
the specified sequence positions. Interactions within and between different
parts of the structure can be individually reweighted, depending on the
desired objective.

Updates to this protocol capture can be found at:
http://kortemmelab.ucsf.edu/data/

Running the Protocol
--------------------
To run this protocol, the backrub app is used to generate an ensemble of
structures. After that, sequence_tolerance is used to sample a large number
of sequence for each member of the ensemble. Finally, an R script is used to
process the results. A python script is included which handles generation of
single ensemble member using backrub and sequence scoring using
sequence_tolerance. It includes a few paths which must be customized to run
correctly on a user's system. Please note that 20 backbones is the minimum
suggested to get acceptable output. The more backbone structures that are
generated, the less prone the results will be to stochastic fluctuations.

Common flags:

    -s
      This flag specifies the starting structure.
    -resfile
      This is used in backrub and sequence_tolerance to specify mutations and 
      control sequence sampling. It is required for sequence_tolerance.
    -score:weights
      This flag is used to specify a weights file that disables environment 
      dependent hydrogen bonds.
    -score:patch
      This flag must be used to reapply the score12 patch to the standard scoring
      function.
    -ex1 -ex2 -extrachi_cutoff
      These flags enable higher resolution rotamer librares for mutation and
      sequence redesign.

Backrub flags:

    -backrub:ntrials
      This flag is used to increase the number of Monte Carlo steps above the
      default of 1000.
    -backrub:minimize_movemap
      If mutations are specified in the resfile, this movemap is used to 
      specify degrees of freedom to be minimized in a three stage process:
      CHI, CHI+BB, CHI+BB+JUMP.
    -in:file:movemap -sm_prob
      Both of these flags are required to enable small phi/psi moves during
      backrub sampling.

Sequence tolerance flags:

    -ms:checkpoint:prefix -ms:checkpoint:interval
      Both of these flags must be specified to get output of the scored sequences.
    -ms:generations -ms:pop_size -ms:pop_from_ss
      These flags affect the genetic algorithm used for sequence sampling.
    -score:ref_offsets
      This flag is used to reweight the reference energies for given residues.
    -seq_tol:fitness_master_weights
      This flag controls the fitness function used for the genetic algorithm.

Example Rosetta command-line:

    rosetta-3.2/rosetta_source/bin/backrub.linuxgccrelease -database rosetta-3.2/rosetta_database -s input_files/2I0L_A_C_V2006/2I0L_A_C_V2006.pdb -ex1 -ex2 -extrachi_cutoff 0 -mute core.io.pdb.file_data -backrub:ntrials 10000 -score:weights input_files/standard_NO_HB_ENV_DEP.wts -score:patch score12
    rosetta-3.2/rosetta_source/bin/sequence_tolerance.linuxgccrelease -database rosetta-3.2/rosetta_database -s 2I0L_A_C_V2006_0001_low.pdb -ex1 -ex2 -extrachi_cutoff 0 -score:ref_offsets HIS 1.2 -seq_tol:fitness_master_weights 1 1 1 2 -ms:generations 5 -ms:pop_size 2000 -ms:pop_from_ss 1 -ms:checkpoint:prefix 2I0L_A_C_V2006_0001 -ms:checkpoint:interval 200 -ms:checkpoint:gz -score:weights input_files/standard_NO_HB_ENV_DEP.wts -out:prefix 2I0L_A_C_V2006_0001 -score:patch score12 -resfile input_files/2I0L_A_C_V2006/2I0L_A_C_V2006_seqtol.resfile

Using the seqtol_resfile.py python script
-----------------------------------------

The seqtol_resfile.py takes as input a PDB file and generates a resfile for
use with the sequence_tolerance app. It takes at least two other required
arguments. The first is the command used for making residues designable. This
is usually either "ALLAA" for all amino acids, or "PIKAA ..." for a 
restricted set of amino acids. The next arguments are the residues which 
should be designable, with the chain and residue number separated by a
colon.

Example seqtol_resfile.py command-line:

    scripts/seqtol_resfile.py input_files/2I0L_A_C_V2006/2I0L_A_C_V2006.pdb "PIKAA ADEFGHIKLMNPQRSTVWY" B:2002 B:2003 B:2004 B:2005 B:2006

Using the backrub_seqtol.py python script
-----------------------------------------

The backrub_seqtol.py script takes as input a PDB file and other similarly
named configuration files, and produces a single backrub ensemble member
along with approximately 10,000 scored sequences on that member. All of the
input files use a base name derived from removing the ".pdb" extension from
the PDB file. For instance, the base name of 1MGF.pdb would be 1MFG.

If you want to use one PDB file with many different input files you can 
specify a different path from which to get the input files.

Required input files:

    <base name>_seqtol.resfile
      This resfile specifies which sequence positions to sample, along with the
      residue positions that should be repacked.

Optional input files:

    <base name>_backrub.resfile
      This resfile specifies which residues should have flexible side chains
      during the backrub run. By default, all side chains are flexible. This file
      can also define mutations that should be made to the input structure prior
      to the backrub simulation.
    <base name>_minimize.movemap
      This file is passed to the -backrub:minimize_movemap flag (see above).
    <base name>_perturb.movemap
      This file is passed to the -in:file:movemap flag (see above). It also sets
      -sm_prob flag to 0.1.

Example overall command-line:

    scripts/backrub_seqtol.py input_files/2I0L_A_C_V2006/2I0L_A_C_V2006.pdb 1

Post-processing with R:

    $ cd output/2I0L_A_C_V2006
    $ R
    > source("rosetta-3.2/rosetta_source/analysis/apps/sequence_tolerance.R")
    > process_specificity()

Generating Smith & Kortemme PLoS One 2011 figures:

    $ cd output
    $ R
    > rosetta_analysis_dir <- "rosetta-3.2/rosetta_source/analysis"
    > source("../scripts/figures.R")

Versions
--------

This protocol used the Rosetta 3.2 release version

Several other scripting/analysis tools were used:

* Python 2.4.3
* R 2.12.1

References
----------

* Smith, C. A. & Kortemme, T. (2010) Structure-Based Prediction of the Peptide 
  Sequence Space Recognized by Natural and Synthetic PDZ Domains. J Mol Biol 
  402, 460-474.
* Smith, C. A. & Kortemme, T. (2011) Predicting the Tolerated Sequences for 
  Proteins and Protein Interfaces Using RosettaBackrub Flexible Backbone 
  Design. PLoS One 10.1371/journal.pone.0020451
