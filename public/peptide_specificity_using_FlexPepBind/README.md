Prediction of peptide binding using FlexPepBind
----------------------------------------------

KEYWORDS: PEPTIDES DOCKING INTERFACES

Predicting the structure of a peptide-protein interaction (and an interaction in general) does not guarantee that the peptide indeed binds to the receptor.
However, an accurate model of a peptide-protein interaction can be used to evaluate the binding ability of a given peptide sequence to a given receptor, e.g. relative to other peptides. We have found that this works particularly well when the characteristic binding features of an interactions are defined and reinforced during modeling: sequences for which low-energy models can be generated that fullfill these constraints are assumed to be binders.

We have developed FlexPepBind - a protocol that uses the FlexPepDock framework to model the structure of a range of different peptide-receptor complexes and identifies binders / non-binders.

NOTE: FlexPepBind is NOT a general protocol - it needs to be calibrated for each system based on a small set of experimentally characterized substrates/non-substrates. Once calibrated, new substrates can be identified on a large scale.

This demo illustrates how to run the FlexPepBind protocol to predict peptide binding for two systems: (1) substrates that are deacetylated by Histone Deacetylase 8 (HDAC8), and (2) substrates that are farnesylated (a farnesyl moiety is attached to their c-terminus) by Farnesyl Transferase (FTase).

Protocol overview
-----------------
The FlexPepBind protocol predicts peptide binding by modeling different peptide sequences onto a template peptide-receptor complex. Given such a peptide-receptor complex, a set of peptide sequences are threaded onto the template peptide backbone and further optimized using FlexPepDock protocols.

System-specific calibration: Depending on the system, either simple minimization (of the full peptide conformation and rigid body orientation, and the receptor side chains) or extensive optmization using the FlexPepDock refinement protocol might be needed. Also, the optimal measure to rank the different peptides has to be determined (i.e., a score that highlights the interface energy; see below).

In this demo we show how to use the minimization only protocol to predict binders and non-binders in a set of peptide sequences.

Running the FlexPepBind protocol
--------------------------------
The two main stages of FlexPepBind are:

 1. Thread the query peptide sequence onto the template peptide backbone (using Rosetta fixbb design).
 2. Optimize the threaded complex structure. We use FlexPepDock protocols (either FlexPepDock refinement or minimization only) for this stage.

To run on the specific systems provided here, go to the relevant directory, and:

 a. Run the fpbind_run.sh script (located in the scripts directory). This will minimize the peptide-protein complexes after threading the peptide sequence (listed in the input_files/peptide.list file) one by one onto the template.
 b. After the run finished, run fpbind_analysis.sh (located in the scripts directory). It will extract the relevant scores of the minimized structures & save it in the score_analysis/ directory.

Scores currently considered are:
 1. Interface score - the energy contributed at the interface (I_sc)
 2. peptide score - the energy contributed by the internal energy of the peptide + the interface score (pep_sc)
 3. peptide score without reference energy - the same peptide score, but without the amino acid dependent energy terms (Eaa) that were optimized to generate designs with natural amino acid content (pepe_sc_noref; we found that removing this term can significantly improve distinction of substrates from non-substrates in certain systems)
 4. Reweighted score: = sum (total_score + peptide_score + interface score) (reweighted_sc)

You need to change the paths of the Rosetta executables and database directories in the fpbind_run script.

 ROSETTA_BIN="rosetta/main/source/bin"

 ROSETTA_DB="rosetta/main/database/"

Example run:

 $ cd HDAC8/
 $ bash ../scripts/fpbind_run.sh
 $ bash ../scripts/fpbind_analysis.sh

The output is the score with different scoring terms:
For example, in score_analysis/I_sc

 $ paste input_files/peptide.list score_analysis/I_sc

 GYKFGC	-16.757

 GFKWGC	-17.108

 GFKFGC	-16.352

 GMKDGC	-13.870

 GDKDGC	-13.179

 GQKIGC	-13.408

The above lines show the interface scores for 6 peptides (the top 3 are HDAC8 substrates & the bottom 3 are not HDAC8 substrates).


Further information
-------------------
London N., Lamphear C.L., Hougland J.L., Fierke C.A., and Schueler-Furman O. (2011).Identification of a novel class of farnesylation targets by structure-based modeling of binding specificity. PLoS Comput Biol, 7: e1002170.

London N., Gulla S., Keating A.E., and Schueler-Furman O. (2012). In silico and in vitro elucidation of BH3 binding specificity toward Bcl-2. Biochemistry, 51: 5841-50.

Alam N., Zimmerman L., Wolfson N.A., Joseph C.G., Fierke C.A., and Schueler-Furman O. (2016). Structure-Based Identification of HDAC8 Non-histone Substrates. Structure, 24: 458-68.
