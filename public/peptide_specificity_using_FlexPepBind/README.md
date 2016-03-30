Peptide specificty using FlexPepBind
------------------------------------
This demo illustrates how to run FlexPepBind protocol to predict peptide binding specificities of HDAC8 and FTase.

Protocol overview
-----------------
The FlexPepBind protocol predicts peptide binding specificity by modeling different peptide sequences onto a template peptide-receptor complex. Given a peptide-receptor complex a set of peptide sequences are threading onto the template peptide backbone and further optimized using FlexPepDock protcols. Depeneding on the system either simple minimization or extensive optmization using the refinement protocol might be needed. Here in this demo we have shown how to use minimization only protcol to predict binders and non-binder from a set of peptide sequences.

Running the FlexPepBind protocol
--------------------------------
The two main stages of FlexPepBind are:
1. Threading the query peptide sequence onto the template peptide backbone. This is done using Rosetta fixbb design protocol
2. Optimizing the threaded complex. We use FlexPepDock protocols (either FlexPepDock refinement or minimization only) for this stage.

To run on the specific systems here provided go to that directory and then

a. Run fpbind_run.sh script located in the scripts directory. This will minimize the peptide-protein complexes after threading peptide sequence (listed in the input_files/peptide.list file) one by one onto the template.
b. After the ran finished run fpbind_analysis.sh located in the scripts directory. It will extract the relevent scores (I_sc, reweighted_sc, pep_sc, pep_sc_noref) of the minimized structures & save it in score_analysis/ directory.

You need to change the paths of the Rosetta executables and database directories in the fpbind_run script.

 ROSETTA_BIN="rosetta/main/source/bin"
 ROSETTA_DB="rosetta/main/database/"

Ex.
$ cd HDAC8/
$ scripts/fpbind_run.sh
$ scripts/fpbind_analysis.sh

The output is the score with different scoring term
For example in score_analysis/I_sc

$ paste input_files/peptide.list score_analysis/I_sc
GYKFGC	-16.757
GFKWGC	-17.108
GFKFGC	-16.352
GMKDGC	-13.870
GDKDGC	-13.179
GQKIGC	-13.408

The above line shown the scores for 6 peptide (top 3 are HDAC8 substrates & bottom 3 are not HDAC8 substrates)


Further information
-------------------
London N., Lamphear C.L., Hougland J.L., Fierke C.A., and Schueler-Furman O. (2011).Identification of a novel class of farnesylation targets by structure-based modeling of binding
specificity. PLoS Comput Biol, 7: e1002170.

London N., Gulla S., Keating A.E., and Schueler-Furman O. (2012). In silico and in vitro elucidation of BH3 binding specificity toward Bcl-2. Biochemistry, 51: 5841-50.

Alam N., Zimmerman L., Wolfson N.A., Joseph C.G., Fierke C.A., and Schueler-Furman O. (2016). Structure-Based Identification of HDAC8 Non-histone Substrates. Structure, 24: 458-68.
