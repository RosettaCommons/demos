# SWA Protein Main

KEYWORDS: STRUCTURE_PREDICTION LOOPS

## Author
Rhiju Das, rhiju@stanford.edu

# Loop Remodeling by Enumeration: the Core Step of Protein 'Stepwise Assembly'

## Brief Description

Build a loop denovo by enumerating through phi,psi angles, and closing the chain by CCD. Should give a complete enumeration for a loop up to 5 residues in length.

## Abstract

Consistently predicting protein structure at atomic resolution from sequence alone remains an unsolved problem in computational biophysics. Practical challenges involving protein loops arise frequently in ab initio modeling, comparative modeling, and protein design, but even these cases can become intractable as loop lengths exceed 10 residues and if surrounding side-chain conformations are erased. This demo illustrates a novel approach to protein modeling that is more powerful than prior methods that strive for atomic resolution. The central innovation is a ‘stepwise ansatz’ inspired by recent ab initio RNA algorithms, which resolves a conformational sampling bottleneck through residue-by-residue conformer enumeration and dynamic programming.


Reference: R. Das (2013) "Atomic-accuracy prediction of protein loop structures enabled by an RNA-inspired ansatz", under review.
More info: http://arxiv.org/abs/1208.2680

## Example Rosetta Command Line

This rebuilds residues 5-8 on the knottin scaffold 2it7 which has had the loop and all the protein's sidechains removed:

```
swa_protein_main -rebuild  -s1 rosetta_inputs/noloop5-8_2it7_stripsidechain.pdb   -input_res1 1-4 9-28   -sample_res 5 6  -bridge_res 7 8  -cutpoint_closed 7    -superimpose_res 1-4 9-28  -fixed_res 1-4 9-28   -calc_rms_res 5-8  -jump_res 1 28  -ccd_close  -out:file:silent_struct_type binary  -fasta rosetta_inputs/2it7.fasta  -n_sample 18  -nstruct 400  -cluster:radius 0.100  -extrachi_cutoff 0  -ex1  -ex2  -score:weights score12.wts  -pack_weights pack_no_hb_env_dep.wts  -in:detect_disulf false  -add_peptide_plane  -native rosetta_inputs/2it7.pdb  -mute all   -out:file:silent 2it7_rebuild.out -disulfide_file rosetta_inputs/2it7.disulf
```

If you plot column 26 against 1 in 2it7_rebuild.out you should see that the lowest energy models have backbone rmsds of well under 1 Angstrom.

```
%(bin)s/extract_pdbs.%(binext)s -in:file:silent 2it7_rebuild.out   -in:file:silent_struct_type binary  -in:file:tags S_0 -database %(database)s -run:constant_seed -nodelay  2>&1 \
```

**NOTE:** Running 'StepWise Assembly' on a longer loop requires a more complex workflow that carries out buildup of the loop across all possible residue-by-residue build paths. This requires a master python script to setup the job, and another master python script to queue up the resulting computation, which is described as a directed acyclic graph. This full workflow is being presented in a separate demo [swa_protein_long_loop].

## Versions
This should work directly out of trunk for any version of Rosetta after June 2012; however for versions before March 2013, the name of the "swa_protein_main" application was "stepwise_protein_test".


