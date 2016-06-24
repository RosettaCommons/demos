#Relax with all-heavy-atom constraints

KEYWORDS: UTILITIES GENERAL

##Introduction

We looked for a way to simultaneously minimize rosetta energy and keep all heavy atoms in a crystal structure as close as possible to their starting positions. As hard-won experience has shown, simply running relax on a structure will often move the backbone a few Angstroms. The best way we have found to perform the simultaneous optimization is to run relax with constraints always turned on (typically constraints ramp down in the late cycles of a relax run) and to constrain not just backbone but also sidechain atoms. This protocol has been tested on a benchmark set of 51 proteins and found to increase sequence recovery in enzyme design by 5% as compared with design in raw pdb structures. It accomplishes this with only .077 Angstrom RMSD over the set of proteins (C-alpha RMSD) from raw pdb to relaxed-with-csts pdb. A more complete description of the data leading to this protocol is below.


## Protocol

### Prepare structures for relax

The required files are in: `rosetta/rosetta_source/src/apps/public/relax_w_allatom_cst`
Many pdbs have features, such as non-canonical amino acids, which will cause rosetta to fail. This script will be able to process most input pdbs so they are ready for a relax run. The script will "clean" structures to replace non-canonical amino acids with their closest counterparts: 
```
rosetta/rosetta_source/src/apps/public/relax_w_allatom_cst/clean_pdb_keep_ligand.py

python rosetta/rosetta_source/srsrc/apps/public/relax_w_allatom_cst/clean_pdb_keep_ligand.py your_structure_original.pdb -ignorechain
```

### Relax with all-heavy-atom constraints: Short protocol (recommended)

Relax with all-heavy-atom constraints is built into the relax application itself. If this is a new structure you may want to first clean it up using the above script. Relax proceeds as follows:
(Note that these flags will not preserve any ligands in the structure unless the ligand params file is added.)
```
$> cp starting_inputs/1A99_1A99.pdb .

$> $ROSETTA3/bin/relax.default.linuxgccrelease -s 1A99_1A99.pdb @starting_inputs/flags2 > log2.txt
```
(where `$ROSETTA3`=path-to-Rosetta/main/source)

These flags are required: 
```
-relax:constrain_relax_to_start_coords
-relax:coord_constrain_sidechains
-relax:ramp_constraints false
```
The flags2 file includes a set of recommended flags:
```
-ignore_unrecognized_res
-relax:constrain_relax_to_start_coords
-relax:coord_constrain_sidechains
-relax:ramp_constraints false
-ex1
-ex2
-use_input_sc
-correct
-no_his_his_pairE
-no_optH false
-flip_HNQ
```

### Relax with all-heavy-atom constraints: Longer Protocol (not recommended)

In general the short protocol is preferred for most applications, since this version is more complicated and the two give nearly identical results. In this protocol a separate script first generates sidechain atom constraints from an input pdb, then the relax protocol is run with this pre-generated constraitn file. The shorter protocol does this all in one step, and this longer version is largely deprecated. Certain users might prefer this protocol because it allows you to see a list of all constraints, and perhaps to modify constraints using other scripts/data, prior to relax. 

(a) Generate sidechain coordinate constraints on your pdb using sidechain_cst_3.py:
```

$> cp starting_inputs/1A99_1A99.pdb .
$> $ROSETTA3/src/apps/public/relax_w_allatom_cst/sidechain_cst_3.py 1A99_1A99.pdb 0.1 0.5
[output:1A99_1A99_sc.cst]
```

(b) Run relax using these constraints and with a custom relax script to force constraints to stay on during the entire run. 
Note that these flags will not preserve any ligands in the structure unless the ligand params file is added. If you want to keep the ligand, simply copy it over from the input pdb. 
```
$> $ROSETTA3/bin/relax.default.linuxgccrelease -s 1A99_1A99.pdb -constraints:cst_fa_file 1A99_1A99_sc.cst @starting_inputs/flags > log.txt
```
flags file:
```
-constrain_relax_to_start_coords
-relax:script ../../../../rosetta_source/src/apps/public/relax_w_allatom_cst/always_constrained_relax_script
-ignore_unrecognized_res
-preserve_header
-ex1
-ex2
-use_input_sc
-correct
-no_his_his_pairE
-score::hbond_params correct_params
-lj_hbond_hdis 1.75
-lj_hbond_OH_donor_dis 2.6
-linmem_ig 10
-nblist_autoupdate true
-dun08 false
-no_optH false
-flip_HNQ
-chemical:exclude_patches LowerDNA  UpperDNA Cterm_amidation SpecialRotamer
 VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals
pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated
thr_phosphorylated  tyr_phosphorylated tyr_sulfated lys_dimethylated
lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated
cys_acetylated tyr_diiodinated N_acetylated C_methylamidated
MethylatedProteinCterm
```
note: Including extra rotamers is important if your goal is to keep all sidechain atoms tightly constrained; if that is not important for your applications, by all means exclude the ex1 and ex2.

### Relax with all-heavy-atom constraints: Data
To test protocols for relaxation of input pdbs we used the 51 scaffold test set for enzdes. We run design over the input structures and calculate the percent of residues which come back with the native identity -- sequence recovery, only over the designed residues. This is of course an imperfect metric (the original sequence might not be fully optimal), but it allows us to ask how many residues rosetta will correctly choose, assuming that the input structure is already at a minimum in sequence space for the ligand in question. It had already been found that relax alone (with no constraints) will distort most structures, and that those structures will give a much higher sequence recovery in design, but this is a result of the distortion to the input structure.

We decided to test a few protocols to find the one that would best minimize RMSD from the original pdb while maximizing sequence recovery. All calculations are averages over 50 runs of the 51 scaffold set. All RMSDs are to native C-alpha. It is also possible to read in electron density for a pdb and use that as a constraint during relax, but electron density is not uniformly available, and we were looking for a protocol that would work in every case even if the electron density was lacking. We chose the first protocol as the best because it minimized rmsd (and sidechain motions) while maximizing sequence recovery.

Running relax with sidechain coordinate constraints and bb coord constraints: (ex flags and use native; ex flags on in enzdes)

    0.447 sequence recovery (0.077 RMSD) [-557 totalscore]

Running relax with sc-sc distance constraints at 3.5 distance cutoff. Note that this gets similar rmsd minimization but doesn't maximize sequence recovery or minimize score as well as coordinate constraints. We tested this protocol with a variety of sc-sc distance constraint cutoff values and found it to slightly but systematically underperform coordinate constraints:

    0.436 sequence recovery (0.0706 RMSD) (-534 totalscore)

Running elax with backbone constraints only. Note that this does worse in terms of rmsd and that many sidechains are a few angstroms off:

    0.488 sequences recovery (0.098 RMSD) (-633 totalscore)

No relax benchmark and rosetta scoring of native input structures:

    0.40 (0 RMSD by definition) (-194.7 totalscore)

At this point, the astute reader might ask, what score terms became 438 Rosetta Energy Units (REU) better? We ranked the difference in scores over all structures, comparing the all-atom coordinate constraint protocol and the non-relaxed input structure. The biggest difference is fa_dun (-192.4), followed by fa_rep (-108.1), pro_close (-26.4), hbond_sc (-10.8) and omega (-8.5). Many input rotamers are close to but not in a "good" Dunbrack rotamer, and the backbone has to be slightly tweaked in order for that residue to get a good dunbrack score. Also many atoms are slightly too close, and they give the fa_rep contribution.


