This file contains information specific for FTase FlexPepBind. For general information about FlexPepBind, see README file in the upper level directory.
KEYWORDS: PEPTIDES DOCKING INTERFACES
Input files needed to run the FlexPepBind protocol on FTase:
------------------------------------------------------------

1tn6.pdb      : The Ftase-peptide complex structure
template.pdb  : The template structure created from 1tn6
substrates    : list of Ftase peptide substrates
nonsubstrates : list of FTase peptide non-substrates
peptide.list  : list of peptide sequence to the the FlexPepBind protocol on. It contains both the substrates and the non-substrates
resfile       : threading instruction (see more on resfile at https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/d1/d97/resfiles.html)
fixbb_flags   : list of flags for running peptdide threading
peptide.list  : list of peptide sequence on which to run the FlexPepBind protocol. It contains both the substrates and the non-substrates
resfile       : threading instructions (see more on resfile at https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/d1/d97/resfiles.html)
minimization_flags : file containing list of flags for running FlexPepDock minimization
fpp.params    : param file for the small non-peptidic molecule bound to FTase at the active site (see more on how to create param files for non-protein molecules here: https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/df/de9/preparing_ligands.html)
constraints.cst : list of constraints used during the minimization


Generation of the template structure (template.pdb) for FTase FlexPepBind:
--------------------------------------------------------------------------
Starting with the complex structure 1tn6, a template pdb was created following the steps mentioned below:
a. The input pdb was cleaned by removing solvent molecules.
b. A param file for FPP was created.
c. The receptor structure was prepacked to remove internal clashes. (see the demo refinement_of_protein_peptide_complex_using_FlexPepDock/ for more information on how to run prepacking). The FPP molecule was read in using the flag -extra_res_fa FPP.params.
d. The prepacked structure was used for running FlexPepBind.


How to run the FTase FlexPepBind protocol:
------------------------------------------

Run the FTase FlexPepBind protocol as:
$ bash ../scripts/fpbind_run.sh
$ bash ../scripts/fpbind_analysis.sh

The scoring term that worked best for FTase FlexPepBind is pep_sc_noref. The output is in score_analysis/pep_sc_noref

$ paste input_files/peptide.list score_analysis/pep_sc_noref

CSII	-4.9

CLIT	-3.2

CFLS	-2.2

CKKP	 7.6

CTKR	10.7

CSIP	 5.4

The pep_sc_noref scores are listed for 6 peptides (the top 3 are FTase substrates & the bottom 3 are not FTase substrates; it might differ depending on the rosetta version and scoring function used to run the protocol)
