Scoring Tutorial
================

KEYWORDS: CORE_CONCEPTS ANALYSIS UTILITIES GENERAL  

Tutorial by Shourya S. Roy Burman (ssrb@jhu.edu) 
Created 20 June 2016

Updated 29 May 2017 by Vikram K. Mulligan (vmullig@uw.edu) for the ref2015 energy function.

[[_TOC_]]

Summary
-------
Rosetta calculates the energy of a biomolecule using energy functions based on physical and statistical terms. By the end of this tutorial, you should understand:

* What it means to calculate the energy of a biomolecule in Rosetta
* How much a particular type of interaction contributes to the energy/score
* How to interpret and compare the energy scores calculated by Rosetta
* How to prepare a biomolecular structure for scoring
* How to score a biomolecule
* How to change the energy function to other preset functions
* How to customize the terms in an energy function for your purpose
* How to get every residue's contribution to the energy score

Scoring in Rosetta
------------------

In Rosetta, the energy of a biomolecule is calculated by _scoring_ it. Rosetta has an optimized energy function or _score function_ called _ref2015_ for calculating the energy of all atomic interactions in a globular protein made of L-amino acids. There are also several all-atom score functions for specialized applications on other biomolecules, as well as score functions for the reduced [[_centroid representation_../full_atom_vs_centroid/fullatom_centroid.md]]. Additionally, you can create a custom score function to suit your requirements.

Score Function
--------------

Score functions in Rosetta are weighted sums of energy terms, some of which represent physical forces like electrostatics and van der Waals' interactions, while others represent statistical terms like the probability of finding the torsion angles in Ramachandran space. Below is a list of the energy terms used in the _ref2015_ score function:

    fa_atr                 Lennard-Jones attractive between atoms in different residues
    fa_rep                 Lennard-Jones repulsive between atoms in different residues
    fa_sol                 Lazaridis-Karplus solvation energy
    fa_intra_sol_xover4    Intra-residue Lazaridis-Karplus solvation energy
    lk_ball_wtd            Asymmetric solvation energy
    fa_intra_rep           Lennard-Jones repulsive between atoms in the same residue
    fa_elec                Coulombic electrostatic potential with a distance-dependent dielectric   
    pro_close              Proline ring closure energy and energy of psi angle of preceding residue
    hbond_sr_bb            Backbone-backbone hbonds close in primary sequence
    hbond_lr_bb            Backbone-backbone hbonds distant in primary sequence
    hbond_bb_sc            Sidechain-backbone hydrogen bond energy
    hbond_sc               Sidechain-sidechain hydrogen bond energy
    dslf_fa13              Disulfide geometry potential
    rama_prepro            Ramachandran preferences (with separate lookup tables for pre-proline positions and other positions)
    omega                  Omega dihedral in the backbone. A Harmonic constraint on planarity with standard deviation of ~6 deg.
    p_aa_pp                Probability of amino acid, given torsion values for phi and psi
    fa_dun                 Internal energy of sidechain rotamers as derived from Dunbrack's statistics
    yhh_planarity          A special torsional potential to keep the tyrosine hydroxyl in the plane of the aromatic ring
    ref                    Reference energy for each amino acid. Balances internal energy of amino acid terms.  Plays role in design.
    METHOD_WEIGHTS         Not an energy term itself, but the parameters for each amino acid used by the ref energy term. 


Further description of energy terms can be found [here](https://www.rosettacommons.org/docs/latest/rosetta_basics/scoring/score-types).

The weights associated with the _ref2015_ score function are:

```
METHOD_WEIGHTS ref 0.773742 0.443793 -1.63002 -1.96094 0.61937 0.173326 0.388298 1.0806 -0.358574 0.761128 0.249477 -1.19118 -0.250485 -1.51717 -0.32436 0.165383 0.20134 0.979644 1.23413 0.162496 
fa_atr 1
fa_rep 0.55
fa_sol 0.9375
fa_intra_rep 0.005
fa_elec 0.875
pro_close 1.25
hbond_sr_bb 1.17
hbond_lr_bb 1.17
hbond_bb_sc 1.17
hbond_sc 1.1
dslf_fa13 1.25
rama 0.25
omega 0.625
fa_dun 0.7
p_aa_pp 0.4
yhh_planarity 0.625
ref 1
```

The algorithm that Rosetta uses to score a structure can be found [here](https://www.rosettacommons.org/docs/latest/rosetta_basics/scoring/scoring-explained).

Comparing Rosetta Scores to Real-Life Energies
----------------------------------------------

While much of the energy function in Rosetta is physics-based, it also has certain statistical terms to favor structures that look like known protein structures (as nature often conserves protein folds).

>**While a lower scoring structure is more likely to be closer to the native structure, the scores do not have a direct conversion to physical energy units like kcal/mol. Instead we represent them in _Rosetta Energy Units (REU)_.**

Also, as the scores depend on the score function used, it is usually not meaningful to compare structures scored using different score functions.

Navigating to the Demos
-----------------------

The demos are available at `<path_to_Rosetta_directory>/demos/tutorials/scoring`. All demo commands listed in this tutorial should be executed when in this directory. All the demos here use the `linuxgccrelease` binary. You may be required to change it to whatever is appropriate given your operating system and compiler.

Demo
----

###Preparing PDBs for Scoring

To score a biomolecule in Rosetta, we use the `score_jd2` executable. This application, along with most in Rosetta, expects the input PDB to be formatted in a certain manner. A PDB downloaded directly from the Protein Data Bank may or may not work with Rosetta in general, and `score_jd2` in particular. Here's an example where we try to score the PDB 3TDM. When in the right demo directory, run:

    > $ROSETTA3/bin/score_jd2.linuxgccrelease -in:file:s input_files/from_rcsb/3tdm.pdb
   
The application will exit with the following error:

    ERROR: Unrecognized residue: PO4
    
This PDB contains a phosphate ion that Rosetta is unable to process without additional options. To score this PDB, we will add an option `-ignore_unrecognized_res`, which simply ignores the phosphates in the PDB.

    $> $ROSETTA3/bin/score_jd2.linuxgccrelease -in:file:s input_files/from_rcsb/3tdm.pdb -ignore_unrecognized_res

Now the PDB will be scored and the score will be displayed in a file `score.sc` (henceforth called the _score file_) in your current working directory. We will learn how to analyze and interpret this file in the next section. *To proceed on to the next step, remove `score.sc` by typing `rm score.sc`. Otherwise, all energy scores of the structures scored here onwards will be appended to this file.*

<a name="output_struc"></a>If an input PDB does not meet the exact specifications that Rosetta requires, e.g. it has missing heavy atoms or unusual residues that Rosetta recognizes by default (unlike phosphates), Rosetta adds or changes atoms to satisfy the specifications. You can ask it to output the structure it actually scores by including the option `-out:pdb`. In this example, we will score the PDB 1QYS taken directly from the Protein Data Bank:

    $> $ROSETTA3/bin/score_jd2.linuxgccrelease -in:file:s input_files/from_rcsb/1qys.pdb -out:pdb
    
In the log, you will see the following lines:

```
...
core.io.pose_from_sfr.PoseFromSFRBuilder: Reading MSE as MET!
...
core.pack.pack_missing_sidechains: packing residue number 13 because of missing atom number 6 atom name  CG
...
```

The first line indicates that it converts the residue _MSE_, i.e. selenomethionine to _MET_, i.e. regular methionine. The second line tells you that Rosetta found that the C<sub>γ</sub> atom was missing in residue number 13, and built the sidechain for residue number 13. You should see that `score.sc` has been written to and a new file `1qys_0001.pdb` is present in your current working directory.

>**Since Rosetta non-deterministically rebuilds missing sidechains, every run of this example will produce a different result, both in terms of the PDB structure and the score file.**

The PDB file now contains the missing atoms and MET in place of MSE. This is the structure that was actually scored. Also, the score file will show a large positive `total_score` indicating an unfavorable structure. **This does not mean that the structure is unstable, it simply means that Rosetta believes that some minor steric clashes may exist in this PDB.** The score file will have a similar format to the following:

```
SEQUENCE: 
SCORE: total_score dslf_fa13    fa_atr    fa_dun   fa_elec fa_intra_rep       fa_rep       fa_sol hbond_bb_sc hbond_lr_bb    hbond_sc hbond_sr_bb linear_chainbreak             omega overlap_chainbreak            p_aa_pp pro_close      rama       ref yhh_planarity description 
SCORE:     267.496     0.000  -422.275   290.201   -25.824        1.313      238.436      248.433      -1.045     -23.835      -2.245     -22.744             0.000             1.234              0.000             -4.258     0.000     2.749   -12.643         0.000 1qys_0001
```

*To proceed on to the next step, remove `score.sc` by typing `rm score.sc`. Otherwise, all energy scores of the structures scored here onwards will be appended to this file.*

To avoid these issues, it is recommended that you always refine the PDB with the [[relaxi|Relax]] protocol with the same _score function_ that you intend to eventually score with. This relieves clashes and prepares the structure for scoring in Rosetta. More details on how to do this in later tutorials. Let us switch our focus to scoring refined PDBs.

###Basic Scoring

In this section, we are going to score the PDB 1QYS (a refined version is provided in `<path_to_Rosetta_directory>/demos/tutorials/scoring/input_files`). First, we will use the default score function, i.e. _ref2015_. Instead of passing a whitespace separated list of options, we will start using flags files as we start using more options. For now, the only options that we will pass in the flags file is the input PDB and the output score file name:

    -in:file:s input_files/1qys.pdb
    
    -out:file:scorefile output_files/score.sc
    
To score this structure, run:

    $> $ROSETTA3/bin/score_jd2.linuxgccrelease @flag
    
Running this should produce a file called `score.sc` in the directory `output_files`. Compare this to the file `<path_to_Rosetta_directory>/demos/tutorials/scoring/output_files/expected_output/score.sc`. They should be the same (except for the `time` column which depends on your CPU speed).

####Analysis of the Score File

The `score.sc` file should look like:

```
SEQUENCE: 
SCORE: total_score       score dslf_fa13    fa_atr    fa_dun   fa_elec fa_intra_rep       fa_rep       fa_sol hbond_bb_sc hbond_lr_bb    hbond_sc hbond_sr_bb linear_chainbreak             omega overlap_chainbreak            p_aa_pp pro_close      rama       ref      time yhh_planarity description 
SCORE:    -163.023    -163.023     0.000  -423.638   109.662   -46.146        1.040       49.117      241.309      -3.934     -26.998     -11.234     -25.491             0.000             4.211              0.000            -13.603     0.000    -4.905   -12.643     1.000         0.230 1qys_0001
```

The first column called `total_score` represents the total weighted score for the structure 1QYS. Notice how much lower the `total_score` is when compared to the unrefined structure in the previous section. For a refined structure of this size, a score of -100 REU to -300 REU is typical. The lower the score, the more stable the structure is likely to be for a given protein.

>**A rule of thumb: -1 to -3 REU per residue is typical while scoring a refined structure with _ref2015_ score function.**

Other columns represent the individual components which go into the total score. For example, the column `fa_atr` represents the weighted score of the Lennard-Jones attractive potential between atoms in different residues. This breakdown can be helpful in determining which energy terms are contributing more than others, i.e. what kind of interactions occur in the protein. In this particular example, Rosetta predicts that the inter-residue van der Waals' forces have the largest stabilizing contribution (`-423.638`), whereas solvation (represented by `fa_sol`) has the highest destabilizing contribution (`241.309`).

>**A large `fa_rep` weighted score (i.e. much larger than stabilizing effect of `fa_atr`) indicates clashes in the structure.**

When we scored the unrefined structure the `fa_rep` weighted score was `238.436`, whereas in the refined structure it is `49.117`. This reduction of ~200 REU indicates that refinement relieved the minor clashes that existed in the original PDB.

###More Scoring Options

####Changing the Score Function

Now, we will try to score two structures using a score function that was optimized to bind two proteins called _docking_. We also want to rename our score file `score_docking.sc`.

A detailed list of score functions can be found in `<path_to_Rosetta_directory>/main/database/scoring/weights`.

The flags file we will use now is:

```
-in:file:l input_files/pdblist

-score:weights docking

-out:file:scorefile output_files/score_docking.sc
```

On running the command

    $> $ROSETTA3/bin/score_jd2.linuxgccrelease @flag_docking

we should get a score file `score_docking.sc` in the directory `output_files`. Compare this to `<path_to_Rosetta_directory>/demos/tutorials/scoring/output_files/expected_output/score_docking.sc`). They should be the same (except for the `time` column which depends on your CPU speed).

```
SEQUENCE: 
SCORE: total_score       score dslf_ca_dih dslf_cs_ang dslf_ss_dih dslf_ss_dst      fa_atr      fa_dun     fa_elec     fa_pair      fa_rep      fa_sol hbond_bb_sc hbond_lr_bb    hbond_sc hbond_sr_bb linear_chainbreak overlap_chainbreak               time description 
SCORE:     -89.596     -89.596       0.000       0.000       0.000       0.000    -143.190       5.640      -1.371      -2.577       3.929      62.290      -0.824      -5.653      -2.502      -5.338             0.000              0.000              0.000 1qys_0001
SCORE:     -55.214     -55.214       0.000       0.000       0.000       0.000    -119.403       4.578      -1.412      -2.044      10.233      66.103      -3.048      -3.799      -3.360      -3.061             0.000              0.000              0.000 1ubq_0001
```

The first line (after the header lines) gives the score of 1QYS and the second line gives the score of 1UBQ. In this score file, we can see that the `total_score` of 1QYS has changed from the previous run. But since we scored the same PDB file, the predicted experimental stability should be the same.

>**It is not meaningful to compare scores generated by different score functions.**

We can also see that `dslf_fa13` energy term from the previous example is missing and is instead replaced by four terms to describe disulfide geometry. The `fa_atr` weighted score has a smaller absolute stabilizing contribution (although the raw, unweighted `fa_atr` score must remain the same, as it is the same structure). The second smaller protein 1UBQ has a higher total score than 1QYS. This is typical for Rosetta, as it calculates per residue scores. This should not be interpreted as 1UBQ being less stable than 1QYS.

>**There does not exist a good correlation between total score and stability of a structure across different proteins.**

####Patch Files and Changing Term Weights

Now, say we want to modify the weights of some of the terms in the _docking_ score function. We may, for example, wish to downweight a more general score term in favor of a more specific set of them, or to add score terms for [constraints](../Constraints_Tutorial/Constraints.md). There are three ways to do this:

* Make a custom weights file, and pass the path to -score:weights
* Add a patch file to modify existing weights
* Set the weight of specific terms from the command line

In the following example, we demonstrate the last two. We use the patch file `docking.wts_patch` which is:

```
fa_atr *= 0.423
fa_rep *= 0.100
fa_sol *= 0.372
fa_intra_rep *= 0.000
fa_pair *= 0.000
fa_dun *= 0.064
hbond_lr_bb *= 0.245
hbond_sr_bb *= 0.245
hbond_bb_sc *= 0.245
hbond_sc *= 0.245
p_aa_pp *= 0.00
fa_elec = 0.026
dslf_ss_dst *= 1.0
dslf_cs_ang *= 1.0
dslf_ss_dih *= 1.0
dslf_ca_dih *= 1.0
pro_close *= 0.000
```

We will also reset the `fa_atr` weight to `1.0` (originally 0.338 in _docking_, and modified to 0.338 * 0.423 = 0.142974 by _docking.wts_patch_) from the command by line using the following options:


```
-in:file:s input_files/pdblist

-score:weights docking
-score:patch docking
-score:set_weights fa_atr 1.0

-out:file:scorefile output_files/score_docking_patch.sc
```

Now on running

    $> $ROSETTA3/bin/score_jd2.linuxgccrelease @flag_docking_patch

we should get a score file `score_docking_patch.sc` in the directory `output_files`. Compare this to `<path_to_Rosetta_directory>/demos/tutorials/scoring/output_files/expected_output/score_docking_patch.sc`). They should be the same (except for the `time` column which depends on your CPU speed).

```
SEQUENCE: 
SCORE: total_score       score dslf_ca_dih dslf_cs_ang dslf_ss_dih dslf_ss_dst      fa_atr      fa_dun     fa_elec      fa_rep      fa_sol hbond_bb_sc hbond_lr_bb    hbond_sc hbond_sr_bb linear_chainbreak overlap_chainbreak               time description 
SCORE:    -404.591    -404.591       0.000       0.000       0.000       0.000    -423.638       0.361      -1.371       0.393      23.172      -0.202      -1.385      -0.613      -1.308             0.000              0.000              1.000 1qys_0001
SCORE:    -332.020    -332.020       0.000       0.000       0.000       0.000    -353.263       0.293      -1.412       1.023      24.590      -0.747      -0.931      -0.823      -0.750             0.000              0.000              0.000 1ubq_0001
```

Note that we have set the weight of `fa_pair` to zero in the patch file. This eliminates the any contribution from that term, and the score file is thus missing the `fa_pair` column. Also, the weighted scores of `fa_rep`, `fa_dun` and others have changed as defined by the patch file. The increase in weighted score of `fa_atr` is a result of the increased weight we fed in through the command line via the flag file. (This, incidentally, is the same weight as the _ref2015_ score function, and hence `fa_atr` has the same weighted score as in the basic scoring example.)

####Advanced Options

Several other options that you could add to the flags file are given [here](https://www.rosettacommons.org/docs/latest/application_documentation/analysis/score-commands).

For example, the option ```-scorefile_format json``` will write the output score file as a JSON file. JSON is a structured data format supported by most modern scripting languages. (For example [Python](https://docs.python.org/2/library/json.html).) This can be useful when writing post-processing scripts.

###Getting Individual Residue Scores

The `total_score` is essentially a sum of the individual residue scores in Rosetta. Sometimes, it might be more useful to get a detailed analysis of how well each residue is scoring. This is especially useful in understanding which residues have clashes in them.

In this section, we will use the executable `per_residue_energies` with the option `-out:file:silent` to specify the file to which to write the per residue breakdown.

```
-in:file:s input_files/1qys.pdb

-out:file:silent output_files/per_res.sc
```

We will run this on the refined 1QYS PDB using

    $> $ROSETTA3/bin/per_residue_energies.linuxgccrelease @flag_per_residue

This should produce a file called `per_res.sc` in `output_files` directory. Compare this to `<path_to_Rosetta_directory>/demos/tutorials/scoring/output_files/expected_output/per_res.sc`). They should be the same. The file will look like:

```
SCORE:     pose_id               pdb_id     fa_atr     fa_rep     fa_sol    fa_intra_rep    fa_elec    pro_close    hbond_sr_bb    hbond_lr_bb    hbond_bb_sc    hbond_sc    dslf_fa13       rama      omega     fa_dun    p_aa_pp    yhh_planarity        ref      score description
SCORE:  input_files/1qys.pdb         3A     -2.666      0.270      2.416           0.025     -0.269        0.000          0.000          0.000          0.000      -0.324        0.000      0.000      0.010      1.894      0.000            0.000     -1.630     -0.273   residue_1
SCORE:  input_files/1qys.pdb         4A     -5.618      0.237      2.802           0.026     -0.123        0.000          0.000         -0.564          0.000       0.000        0.000     -0.262      0.007      0.814     -0.348            0.000      1.081     -1.949   residue_2
...
```

The first line (after the header line) tells us the score breakdown of residue 3 in chain A, represented here as `3A`. (As indicated by the column heading (`pdb_id`), this is in [PDB numbering](Core_Concepts). PDB residue 3A is residue 1 under Rosetta numbering, as indicated by the "residue_1" in the description column.) The output file can be quite difficult to read in text editors because the column headers are often not directly above the column values. You may want open it in MS Excel, Libre Office Calc or a similar spreadsheet program to read them better. 

> **Be sure to use the same score function in `per_residue_energies` as you did in `score_jd2` to have comparable results.**

###Getting Individual Residue Score Breakdowns

If we want to further decompose these scores into inherent residue energies (one-body scores) and residue interaction scores, we will use the `residue_energy_breakdown` with the following options file:

```
-in:file:s input_files/1qys.pdb

-out:file:silent output_files/energy_breakdown.sc
```

We will run this on the refined 1QYS PDB using

    $> $ROSETTA3/bin/residue_energy_breakdown.linuxgccrelease @flag_residue_energy_breakdown

This should produce a file called `energy_breakdown.sc` in `output_files` directory. Compare this to `<path_to_Rosetta_directory>/demos/tutorials/scoring/output_files/expected_output/energy_breakdown.sc`). They should be the same. The file will look like:

```
SCORE:     pose_id                resi1     pdbid1    restype1      resi2     pdbid2    restype2     fa_atr     fa_rep     fa_sol    fa_intra_rep    fa_elec    pro_close    hbond_sr_bb    hbond_lr_bb    hbond_bb_sc    hbond_sc    dslf_fa13       rama      omega     fa_dun    p_aa_pp    yhh_planarity        ref      total                    description
SCORE:  input_files/1qys.pdb          1         3A         ASP         --         --     onebody      0.000      0.000      0.000           0.025      0.000        0.000          0.000          0.000          0.000       0.000        0.000      0.000
...
SCORE:  input_files/1qys.pdb          1         3A         ASP          2         4A         ILE     -1.518      0.072      1.027           0.000      0.721        0.000          0.000          0.000          0.000       0.000        0.000      0.000      0.000      0.000      0.000            0.000      0.000      0.301 input_files/1qys.pdb_1_2
...
```

The first line (after the header line) tells us the score breakdown of the internal (onebody) energies of the first residue (PDB numbering 3A). Further down the file, a line indicates the interaction energy of residues 3A and 4A (shown above). 

> **Be sure to use the same score function in `residue_energy_breakdown` as you did in `score_jd2` to have comparable results.**

###Advanced Scoring Protocols

Some proteins require additional files for scoring. For example, a membrane file requires a definition of the membrane spanning region, while a symmetric protein (using the symmetry framework) requires a definition of the symmetry. Such protocols are beyond the scope of this tutorial. Here are links to the documentation for these methods:

* [Scoring membrane proteins](https://www.rosettacommons.org/docs/latest/application_documentation/membrane_proteins/RosettaMP-App-MPScoring)

* [Symmetric scoring](https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/symmetry#implementation-of-symmetry_symmetric-scoring)

Tips
----

* While the scoring step is deterministic and should provide the same score for given score function and input structure, you may not get the same score if you run `score_jd2` on the same PDB. For example, there might be missing sidechain atoms which Rosetta tries to model (non-deterministically), thus, producing different scores on different runs.
* Do not compare scores produced by different score functions. They could mean very different things.
* Structures should be [relaxed](Relax) into the same score function prior to comparison.

References
----------

A good summary of scoring in Rosetta may be found here:

[Alford et al. (2017) J. Chem. Theory Comput. In press.](https://www.ncbi.nlm.nih.gov/pubmed/28430426)

The default score function in Rosetta, _ref2015_ and its corrections were tested in the paper:

[Park et al. (2016). J. Chem. Theory Comput. 12(12):6201-12](https://www.ncbi.nlm.nih.gov/pubmed/27766851)

The previous default score function, _talaris2014_ was tested in the paper:

[O’Meara et al., J. Chem. Theory Comput. 2015](https://dx.doi.org/10.1021/ct500864r)  

A full description of the changes in its predecessor energy function, _talaris2013_ introduces can be found [here](https://www.rosettacommons.org/node/3508#comment-6946).
