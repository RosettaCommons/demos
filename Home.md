Rosetta Tutorials, Demos, and Protocol Captures
===============================================

[[_TOC_]]

Obtaining tutorial materials
----------------------------

The demos, tutorials, protocol captures, and all example inputs are provided with the full Rosetta distribution, under the demos/ directory. Rosetta is available for license (which is free of charge to academic users) at <https://www.rosettacommons.org/software>. 

<!--- BEGIN_INTERNAL -->
For RosettaCommons users, the demos repository should be automatically downloaded by the get_rosetta.sh download script. Alternatively, RosettaCommons users can download the demos repository from GitHub. e.g.

    git clone git@github.com:RosettaCommons/demos.git 

<!--- END_INTERNAL -->

Tutorials
---------

These are introductory tutorials intended as a gentle introduction to Rosetta concepts, and using common functionality of Rosetta. For additional examples and information on using Rosetta, see the demos (below) or the [Rosetta documentation](https://www.rosettacommons.org/docs/latest/)

Full input files for the tutorials are located in the `demos/tutorials/` directory of the Rosetta distribution. 

### Introduction to Rosetta

1. [[How To Read These Tutorials|Tutorial_Setup]]
2. [[Installing and Building Rosetta|install_build]]
3. [[Working With Rosetta|working_with_rosetta]]
4. [[Controlling Input and Output|input_and_output]]
5. [[Core Rosetta Concepts|Core_Concepts]]
6. [[Working with Non-protein Residues|prepare_ligand_tutorial]] 
7. [[Scoring|scoring]]: Calculating the Energy of a Structure
8. [[Full-Atom vs. Centroid Representations|fullatom_centroid]]
9. [[The Packer|Optimizing_Sidechains_The_Packer]]: Optimizing Sidechains
10. [[Minimization]]: Finding Deeper Energy Wells
11. [[Relax|Relax]]: Refining Structures
12. [[Constraints]]: Biasing Towards a Structure
13. [[Analyzing Rosetta Output|Analysis]]
14. [[The Fold Tree|fold_tree]]: Propagating Changes in the Structure
15. [[Symmetry|Symmetry]]: Modeling Symmetric Proteins
16. [[Scripting with RosettaScripts|scripting_with_rosettascripts]]
17. [[Advanced Scripting with RosettaScripts|advanced_scripting_with_rosettascripts]]
18. [[Commonly Used Options|commonly_used_options]]
19. [[Tips]]

### Commonly Used Rosetta Protocols

* [[de novo (ab initio) Structure Prediction|Denovo_structure_prediction]]
    * [[Advanced de novo Structure Prediction|folding_tutorial]]
* [[Comparative Modeling|rosetta_cm_tutorial]]: Modeling based on Homologs
* Generalized Kinematic Closure (GenKIC): Rapid, versatile loop closure without fragments
    * [[GenKIC Tutorial 1|generalized_kinematic_closure_1]]: Building and closing new loops
    * [[GenKIC Tutorial 2|generalized_kinematic_closure_2]]: Perturbing existing loops
    * [[GenKIC Tutorial 3|generalized_kinematic_closure_3]]: Using pre-selection movers within GenKIC
    * [[GenKIC Tutorial 4|generalized_kinematic_closure_4]]: Closing through disulfides
* [[Loop Modeling and Rebuilding|loop_modeling]]: Modeling Short Fragments
* [[Protein Design|protein_design_tutorial]]
* [[Protein-Protein Docking|Protein-Protein-Docking]]: Modeling Protein-Protein Binding
    * [[Advanced Protein-Protein Docking|advanced_protein-protein_docking_tutorial]]
* [[Protein-Ligand Docking|ligand_docking_tutorial]]: Modeling Protein-Ligand Binding

Demos
-----

Demos are designed to guide users through sample procedures in computational modeling from the point of view of solving a specific problem. 

Full input files for the demos are located in the `demos/public/` directory of the Rosetta distribution.

* [[Demos listed by category|demos-by-category#demos]]
* [[Demos listed by keyword|tag-search]]

Protocol Captures
-----------------

Many papers using Rosetta are accompanied by a protocol capture - an example of how to use the protocol discussed in the paper. The protocol captures below aren't meant to show the best way to solve problems in the current version of Rosetta, instead they are meant to show published solutions to problems that were addressed by members of the Rosetta community. The purpose of these protocol captures is both to serve as a historical record and to assist those trying to reproduce past results. See the demos (above) for updated versions of most protocol captures.

Full input files for the protocol captures are located in the demos/protocol_capture/ directory of the Rosetta distribution.

* [[Protocol captures listed by category|demos-by-category#protocol_captures]]

<!--- BEGIN_INTERNAL --->

Adding new demos
----------------

See the [[How To Write Demos and Tutorials|How_To_Write_Demos_and_Tutorials]] page for details about
writing new demos.

Demos Under Development
-----------------------

If you want to prevent a demo from being published to the static wiki until you've finished developing it or written a paper on it or something, put it in `demos/pilot`.

* [[Devel demos listed by category|demos-by-category#under_development]]

<!--- END_INTERNAL --->
