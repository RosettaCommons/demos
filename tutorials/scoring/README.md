Scoring Tutorial
================
Tutorial by Shourya S. Roy Burman   
Created 20 June 2016

[[_TOC_]]

Summary
-------
Rosetta calculates the energy of a biomolecule using energy functions, some of which are general, while others have been optimized for a particular task or biomolecule. By the end of this tutorial, you should understand:

* What it means to calculate the energy of a biomolecule in Rosetta
* How to interpret and compare the energy scores calculated by Rosetta
* How to score a biomolecule
* How to change the energy function to other preset functions
* How to customize the terms in an energy function for your purpose

Scoring in Rosetta
------------------

In Rosetta, the energy of a biomolecule is calculated by _scoring_ it. Rosetta has an optimized energy function or _score function_ called _talaris2014_ for calculating the energy of all atomic interactions in a globular protein made of L-amino acids. There are also several all-atom score functions for specialized applications on other biomolecules as well as score functions for the reduced _centroid_ represenation. Additionally, you can create a custom score function to suit your requirements.

Score Function
--------------
Score functions in Rosetta are weighted sums of energy terms, some of which represent a physical force like electrostatics and van der Waals', while others represent a statisticial term like the probability of finding the torsion angles in Ramamchandran space. Below is a list of the energy terms used in the _talaris2014_ score function:

    fa_atr                 Lennard-Jones attractive between atoms in different residues
    fa_rep                 Lennard-Jones repulsive between atoms in different residues
    fa_sol                 Lazaridis-Karplus solvation energy
    fa_intra_rep           Lennard-Jones repulsive between atoms in the same residue
    fa_elec                Coulombic electrostatic potential with a distance-dependent dielectric   
    pro_close              Proline ring closure energy and energy of psi angle of preceding residue
    hbond_sr_bb            Backbone-backbone hbonds close in primary sequence
    hbond_lr_bb            Backbone-backbone hbonds distant in primary sequence
    hbond_bb_sc            Sidechain-backbone hydrogen bond energy
    hbond_sc               Sidechain-sidechain hydrogen bond energy
    dslf_fa13              Disulfide geometry potential
    rama                   Ramachandran preferences
    omega                  Omega dihedral in the backbone. A Harmonic constraint on planarity with standard deviation of ~6 deg.
    fa_dun                 Internal energy of sidechain rotamers as derived from Dunbrack's statistics
    p_aa_pp                Probability of amino acid at Φ/Ψ
    ref                    Reference energy for each amino acid. Balances internal energy of amino acid terms.  Plays role in design.
    METHOD_WEIGHTS         Not an energy term itself, but the parameters for each amino acid used by the ref energy term. 


Further description of energy terms can be found [here](https://www.rosettacommons.org/docs/latest/rosetta_basics/scoring/score-types).

The weights associated with the _talaris2014_ score function are:

```html
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

>While a lower scoring structure is more likely to be closer to the native structure, the scores do not have a direct convertion to physical energy units like (kcal/mol). Instead we represent them in _Rosetta Energy Units (REU)_.

Also, as the scores depend on the score function used, it is not fair to compare structures scored using different score functions.

Navigating to the Demos
-----------------------
The demos are available at: `Rosetta/demos/tutorials/scoring`. All demo commands listed in this tutorial should be executed when in this directory. All the demos here use the `linuxgccrelease` binary. You may be required to change it to whatever is appropriate given your operating system and compiler.

Demo
----

###Basic Scoring
In this tutorial, we are going to score the PDB 1QYS (a refined version is provided in `<path_to_Rosetta_directory>/demos/tutorials/scoring/input_files`). First, we will use the default score function, i.e. talaris2014.

    $> <path_to_Rosetta_directory>/main/source/bin/score_jd2.linuxgccrelease @flag
    
The only option that we will pass in the flags file is the input PDB:

    -in:file:s input_files/1qys.pdb
    
    -out:file:scorefile output_files/score.sc
    
Running this should produce a file called `score.sc` in the directory `output_files`. Compare this to the file `<path_to_Rosetta_directory>/demos/tutorials/scoring/output_files/expected_output/score.sc`. They should be the same.

####Analysis of the Score File
The `score.sc` file should look like:

```html
SEQUENCE: 
SCORE: total_score       score dslf_fa13    fa_atr    fa_dun   fa_elec fa_intra_rep       fa_rep       fa_sol hbond_bb_sc hbond_lr_bb    hbond_sc hbond_sr_bb linear_chainbreak             omega overlap_chainbreak            p_aa_pp pro_close      rama       ref      time yhh_planarity description 
SCORE:    -163.023    -163.023     0.000  -423.638   109.662   -46.146        1.040       49.117      241.309      -3.934     -26.998     -11.234     -25.491             0.000             4.211              0.000            -13.603     0.000    -4.905   -12.643     1.000         0.230 1qys_0001
```

The first column called `total_score` represents the total weighted score for the structure 1QYS. For a refined structure of this size, a score of -100 REU to -300 REU is typical.

>A rule of thumb: -2 REU per resdiue is typical.

The column `fa_atr` represents the weighted score of the Lennard-Jones attractive potential between atoms in different residues, and so on. This breakdown can be helpful in determining which energy terms are contributing more than others, i.e. what kind of interactions occur in the protein.

###More Scoring Options
####Changing the score function
Now, we will try to score two structures using an older score function called _score12_. We also want to rename our score file _score_score12.sc.

A detailed list of score functions can be found in `<path_to_Rosetta_directory>/main/database/scoring/weights`.

The flags file we will use now is:

```html
-in:file:l input_files/pdblist

-score:weights score12

-out:file:scorefile output_files/score_score12.sc
```

On running the command

    $> <path_to_Rosetta_directory>/main/source/bin/score_jd2.linuxgccrelease @flag_score12

We should get a score file `score_score12.sc` in the directory `output_files`. Compare this to `<path_to_Rosetta_directory>/demos/tutorials/scoring/output_files/expected_output/score_score12.sc`). They should be the same.

```html
SEQUENCE: 
SCORE: total_score       score dslf_ca_dih dslf_cs_ang dslf_ss_dih dslf_ss_dst      fa_atr      fa_dun fa_intra_rep      fa_pair       fa_rep       fa_sol hbond_bb_sc hbond_lr_bb    hbond_sc hbond_sr_bb linear_chainbreak             omega overlap_chainbreak            p_aa_pp pro_close      rama       ref      time description 
SCORE:    -142.718    -142.718       0.000       0.000       0.000       0.000    -338.910      87.730        0.832       -7.700       39.293      167.307      -3.934     -26.998     -11.234     -12.745             0.000             3.368              0.000            -10.882     0.000    -3.924   -24.920     1.000 1qys_0001
SCORE:     -28.876     -28.876       0.000       0.000       0.000       0.000    -282.611      71.216        0.541       -6.108      102.326      177.549     -14.555     -18.144     -15.085      -7.308             0.000             3.964              0.000             -9.016     0.051    -7.703   -23.990     1.000 1ubq_0001
```

In this score file, we can see that the `total_score` of 1QYS has changed from the previous run. We can also see that `dslf_fa13` energy term from the previous example is missing and is instead replaced by four terms to describe disulphide geometry.

####Patch files and changing term weights
Now, say we want to modify the weights of some of the terms in the score function, score12. There are two ways to do this:
* Adding a patch file providing a list of weights
* Seting the weight of specific terms from the command line

In the following example, we do both. We use the patch file `score12_w_corrections.wts_patch` which is:

```html
hbond_sr_bb *= 0.5
p_aa_pp *= 0.5
rama = 0.2
omega = 0.5
ch_bond_bb_bb = 0.5
fa_cust_pair_dist = 1.0
```
We then set the `fa_atr` weight to `1.0` (originally `0.8` in score12) from the command line using the following flags:


```html
-in:file:s input_files/pdblist

-score:weights score12
-score:patch score12_w_corrections
-score:set_weights fa_atr 1

-out:file:scorefile output_files/score_score12_patch.sc
```

Now on running

    $> <path_to_Rosetta_directory>/main/source/bin/score_jd2.linuxgccrelease @flag_score12_patch

We should get a score file `score_score12_patch.sc` in the directory `output_files`. Compare this to `<path_to_Rosetta_directory>/demos/tutorials/scoring/output_files/expected_output/score_score12_patch.sc`). They should be the same.

```html
SEQUENCE: 
SCORE: total_score       score ch_bond_bb_bb dslf_ca_dih dslf_cs_ang dslf_ss_dih dslf_ss_dst      fa_atr fa_cust_pair_dist            fa_dun fa_intra_rep      fa_pair       fa_rep       fa_sol hbond_bb_sc hbond_lr_bb    hbond_sc hbond_sr_bb linear_chainbreak             omega overlap_chainbreak            p_aa_pp pro_close      rama       ref      time description 
SCORE:    -222.439    -222.439        -6.807       0.000       0.000       0.000       0.000    -423.638             0.000            87.730        0.832       -7.700       39.293      167.307      -3.934     -26.998     -11.234      -6.373             0.000             3.368              0.000             -5.441     0.000    -3.924   -24.920     1.000 1qys_0001
SCORE:     -95.233     -95.233        -3.867       0.000       0.000       0.000       0.000    -353.263             0.000            71.216        0.541       -6.108      102.326      177.549     -14.555     -18.144     -15.085      -3.654             0.000             3.964              0.000             -4.508     0.051    -7.703   -23.990     0.000 1ubq_0001
```
Note that we have new energy terms `ch_bond_bb_bb` and `fa_cust_pair_dist` which were defined in the patch file. Also, the weighted scores of `hbond_sr_bb`, `rama` and others have changed as defined by the patch file. The decrease in weighted score of `fa_atr` is a result of the increased weight we fed in through the command line via the flag file. (This incidentally is the same weight as the _talaris2014_ score function, and hence `fa_atr` has the same weighted score as in the basic scoring example.)

####Advanced options
Several other optons that you could add to the flag file are given [here](https://www.rosettacommons.org/docs/latest/application_documentation/analysis/score-commands).

Tips
----
* While the scoring step is deterministic and should provide the same score for given score function and input structure, you may not get the same score if you run `score_jd2` on PDBs which have not been prepared well. For example, there might be missing sidechain atoms which Rosetta tries to model (non-deterministically), thus, producing different scores on different runs.
* Do not compare scores produced by different score functions. They could mean very different things.
