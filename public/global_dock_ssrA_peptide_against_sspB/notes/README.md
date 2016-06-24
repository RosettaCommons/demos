Author: Barak Raveh
KEYWORDS: PEPTIDES DOCKING
Overview:
=========
This demo demonstrates ab-initio folding and docking of the peptide ssrA (chain D) to the receptor protein sspB (chain B) at a specified binding sites. The user supplies Rosetta with a model of the sspB receptor protein (bound model, for the purpose of this demo), the sequence of the ssrA peptide, and a list of atoms on the surface of sspB that are known to bind the peptide.

Note on identifying the binding site:
================================================
At this stage, Rosetta FlexPepDock allows ab-initio folding and docking of a peptide only at a specifie\
d binding site (completely blind global docking is under development in the Furman lab). That said, for\
tunately peptides have a nice property - many of them tend to bind to large grooves and pockets over th\
e receptor protein surface (London et al., Structure 2010 - The Structural Basis of Peptide-Protein Bin\
ding Strategies). While this does not hold true in 100% of the cases, this is certainly the case for th\
e provided example. Therefore, this demo provides a good example for ab-initio folding and docking of a\
 peptide at a binding site that is specified as a list of atoms or residues over the surface of the rec\
eptor protein. In this case, the binding site residues (specified in <input_files/site_constraints.cst)\
 is based on the Pocket-Finder server (http://www.modelling.leeds.ac.uk/pocketfinder/), which uses the \
Ligsite pocket detection program (Hendlich et al., 1997). Note that other pocket detection prorgrams mi\
ght also work for the same purpose. This would not work in all cases, but here the binding pocket is well defined.

General Sceheme:
================
In order to identify the peptide binding site at the largest pocket, and fold-and-dock the peptide at the binding site, we worked as follows:
(1) We used the Pocket-Finder server (http://www.modelling.leeds.ac.uk/cgi-bin/pocketfinder/pfmage.cgi) on the receptor chain (chain B), and identified the largest pocket. In principle, any other pocket detector can be used for this task. The resulting residues list copied from the server is found in <input_files/pocket_finder_site.txt>.
(2) Based on the site identified in step 1, we created a constraint file that biases FlexPepDock ab-initio to sample peptide conformations at the vicinity of the identified pocket, using soft constraints that can be violated quite generously, to prevent over-bias in case of noise in the input. The constraints for side-chain atoms were repklaced with constraints on the centroid atom, since these constraints are used in the centroid low-resolution step of FlexPepDock. The resulting constrains file is found in <input_files/site_constraints.cst>.
(3) We ran FlexPepDock ab-initio + FlexPepDock Refinement to create the high-resolution model of the interaction. In the starting model, we positioned the peptide in an arbitrary initial orientation relative to the receptor protein, and set its backbone conformation to an ideal extended backbone conformation. The constraint file was used in the initial low-resolution stage of FlexPepDock ab-initio to force the peptide to contact the receptor at the vicinity of the specified binding pocket. 

Details:
========
