Scoring Tutorial
================

In Rosetta, the energy of a biomolecule is calculated by _scoring_ it. Rosetta has an optimized energy function or _score function_ called _talaris2014_ for calculating the energy of all atomic interactions in a globular protein made of L-amino acids. There are also several all-atom score functions for specialized applications on other biomolecules as well as score functions for the reduced _centroid_ represenation. Additionally, you can create a custom score function to suit your requirements.

Score Function
--------------
Score functions in Rosetta are weighted sums of energy terms, some of which represent a physical force like electrostatics and van der Waals', while others represent a statisticial term like the probability of finding the torsion angles in Ramamchandran space. Below is a list of the energy terms used in the talaris2014 score function:


    fa_atr                 Lennard-Jones attractive between atoms in     different residues
    fa_rep                 Lennard-Jones repulsive between atoms in         different residues
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

The algorithm that Rosetta uses to score a structure can be found [here](https://www.rosettacommons.org/docs/latest/rosetta_basics/scoring/scoring-explained).


Demo
----
In this tutorial, we are going to score the PDB 1QYS (a refined version is provided in `Rosetta/demos/tutorials/scoring/input_files`). First, we will use the default score function, i.e. talaris2014.

    Rosetta/main/source/bin/scoring_jd2.linuxgccrelease @flag
    
This should produce a file called `score.sc` in the current working directory. Compare this to the file 'Rosetta/demos/tutorials/scoring/output_files/score.sc'. They should be the same.

Analysis
--------
The `score.sc` file should look like:

```html
SEQUENCE: 
SCORE: total_score       score dslf_fa13    fa_atr    fa_dun   fa_elec fa_intra_rep       fa_rep       fa_sol hbond_bb_sc hbond_lr_bb    hbond_sc hbond_sr_bb linear_chainbreak             omega overlap_chainbreak            p_aa_pp pro_close      rama       ref      time yhh_planarity description 
SCORE:    -163.023    -163.023     0.000  -423.638   109.662   -46.146        1.040       49.117      241.309      -3.934     -26.998     -11.234     -25.491             0.000             4.211              0.000            -13.603     0.000    -4.905   -12.643     1.000         0.230 1qys_0001
```

Flags
-----
    -in:file:fasta ./input_files/1elwA.fasta
    -in:file:frag3 ./input_files/aa1elwA03_05.200_v1_3
    -in:file:frag9 ./input_files/aa1elwA09_05.200_v1_3
    -in:file:native ./input_files/1elw.pdb


    -abinitio:relax
    -nstruct 1
    -out:pdb

    -use_filters true
    -psipred_ss2 ./input_files/1elwA.psipred_ss2
    -abinitio::increase_cycles 10
    -abinitio::rg_reweight 0.5
    -abinitio::rsd_wt_helix 0.5
    -abinitio::rsd_wt_loop 0.5
    -relax::fast

Analyze Output
--------------
The output_files directory contains example output.

In `score.fsc` get a score and RMS for each model.
Typical analysis makes a scatter plot of these with RMS on the x-axis and score on the y-axis.
Look for a "funnel" to low energies and RMS in a successful ab initio prediction.
A failed prediction will not have low RMS/energy structures.
For the following arguments for full production run (using the minirosetta compile):

    -abinitio::fastrelax
    -abinitio::increase_cycles 10
    -abinitio::rg_reweight 0.5
    -abinitio::rsd_wt_helix 0.5
    -abinitio::rsd_wt_loop 0.5
    -abinitio::use_filters false
    -database minirosetta_database
    -ex1
    -ex2aro
    -frag3 aa0000103_05.200_v1_3
    -frag9 aa0000109_05.200_v1_3
    -in:file:boinc_wu_zip ploop23_control_fold_data.zip
    -in:file:native 00001.pdb
    -mute all
    -mute all
    -out:file:silent default.out
    -relax::default_repeats 15
    -silent_gz

    resultfiles = default.out.gz queue = 3000

And run a relax:

    -database
    -ex1
    -ex2aro
    -frag3 aa0000103_05.200_v1_3
    -frag9 aa0000109_05.200_v1_3
    -in:file:boinc_wu_zip ploop23_control_fold_data.zip 
    -in:file:fullatom
    -in:file:native 00001.pdb
    -in:file:s 00001.pdb
     minirosetta_database
    -mute all
    -out:file:silent default.out
    -relax::default_repeats 15
    -run:protocol relax
    -silent_gz

    resultfiles = default.out.gz
  
Tips
----
* While the scoring step is deterministic and should provide the same score for given score function and input structure, you may not get the same score if you run `score_jd2` on PDBs which have not been prepared well. For example, there might be missing sidechain atoms which Rosetta tries to model (non-deterministically), thus, producing different scores on different runs.
