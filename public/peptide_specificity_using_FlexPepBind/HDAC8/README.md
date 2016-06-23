This file contains information specific for HDAC8 FlexPepBind. For general information about FlexPepBind, see README file in the upper level directory.
KEYWORDS: PEPTIDES DOCKING INTERFACES
Input files needed to run the FlexPepBind protocol on HDAC8:
------------------------------------------------------------
template.pdb  : The template structure created from 2v5w
substrates    : list of HDAC8 peptide substrates
nonsubstrates : list of HDAC8 peptide non-substrates
peptide.list  : list of peptide sequence to the the FlexPepBind protocol on. It contains both the substrates and the non-substrates
activity      : percent deacetylation of the peptides listed in peptide.list file
resfile       : threading instruction (see more on resfile at https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/d1/d97/resfiles.html)
fixbb_flags   : file containing list of flags for peptide threading
minimization_flags : file containing list of flags for running FlexPepDock minimization
constraints.cst : list of constraints used during the minimization

Generation of the template structure (template.pdb) for HDAC8 FlexPepBind:
--------------------------------------------------------------------------
The solved structure of a complex of HDAC8 bound to a substrate (PDB id 2v5w) contains only a short peptide-coumarine ligand. Therefore, we needed to generate a structure to serve as template for hexamer-HDAC8 interactions. It was generated following the steps mentioned below:
a. The input pdb was cleaned by removing solvent molecules.
b. Only a monomer bound to a peptide was copied from the solved dimeric structure.
c. The receptor structure was prepacked to remove internal clashes (see the demo refinement_of_protein_peptide_complex_using_FlexPepDock/ for more information on how to run prepacking).
d. The peptide substrates GYKacFGC was threaded onto the template peptide backbone (missing residues were added in extended conformation), and optimize extensively using the FlexPepDock refinement protocol (see the demo refinement_of_protein_peptide_complex_using_FlexPepDock/ for more information on how to run FlexPepDock refinement), and the top-scoring optmized structure (according to interface score; I_sc) was used as the template.

The underlying assumption was that a structural model of a strong substrate (GYKacFGC) would be a good representative of the preferred binding conformation of HDAC8 substrates. Therefore, we used it as a template to evaluate binding ability of other peptides.

How to run the HDAC8 FlexPepBind protocol:
------------------------------------------
$ bash ../scripts/fpbind_run.sh
$ bash ../scripts/fpbind_analysis.sh

The output contains scores for each peptide, according to different scoring terms. For HDAC8 FlexPepBind, we found that the interface score (I_sc) provides the best ranking.

$ paste input_files/peptide.list score_analysis/I_sc
GYKFGC	-16.757
GFKWGC	-17.108
GFKFGC	-16.352
GMKDGC	-13.870
GDKDGC	-13.179
GQKIGC	-13.408

The Interface scores are given for 6 peptides (the top 3 are HDAC8 substrates & the bottom 3 are not HDAC8 substrates; the score might vary depending on the rosetta version and socirng function being used).
