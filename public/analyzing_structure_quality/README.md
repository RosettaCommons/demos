Analyzing Structure Quality
===========================
One of the simplest initial steps to evaluate structural quality is to score the structure using the Rosetta score function and to evaluate the results.

Scoring structures
------------------
To score structures, use the score application. This application can accept any input structure that Rosetta recognizes.

    /path/to/rosetta/bin/score.macosgccrelease -database /path/to/minirosetta_database/ -s ../starting_files/1ubq.pdb.gz -in:file:fullatom -ignore_unrecognized_res -out:output

Description of command-line options
-----------------------------------
* `-database`: Location of the Rosetta database on your system
* `-s`: One or more input strucutre to score, in PDB format.  To use a silent file as input, use `-in:file:silent` instead of `-s` (optionally with `-in:file:tags`).
* `-in:file:fullatom`: Score the structure using the full atom energy function.
* `-out:file:score`: Name of the summary scorefile.
* `-out:output`: Force production of the scored PDBs.

Description of score terms
--------------------------
Rosetta scoring values are reported in REU (Rosetta Energy Units). One study reported a conversion factor of 0.58 from REU to kcal/mol (Kellogg, Leaver-Fay, Baker; Proteins 2011). This conversion may not apply for other applications. Lower (more negative) values indicate a more favorable conformation.

* `fa_atr`:         Lennard-Jones attractive
* `fa_rep`:         Lennard-Jones repulsive
* `fa_sol`:         Lazaridis-Karplus solvation energy
* `fa_intra_rep`:   Lennard-Jones repulsive between atoms in the same residue
* `fa_pair`:        statistics based pair term, favors salt bridges
* `fa_plane`:       pi-pi interaction between aromatic groups, by default = 0
* `fa_dun`:         internal energy of sidechain rotamers as derived from Dunbrack's statistics
* `ref`:            reference energy for each amino acid
* `hbond_lr_bb`:    backbone-backbone hbonds distant in primary sequence
* `hbond_sr_bb`:    backbone-backbone hbonds close in primary sequence
* `hbond_bb_sc`:    sidechain-backbone hydrogen bond energy
* `hbond_sc`:       sidechain-sidechain hydrogen bond energy
* `p_aa_pp`:        Probability of amino acid at phipsi
* `dslf_ss_dst`:    distance score in current disulfide
* `dslf_cs_ang`:    csangles score in current disulfide
* `dslf_ss_dih`:    dihedral score in current disulfide
* `dslf_ca_dih`:    ca dihedral score in current disulfide
* `pro_close`:      proline ring closure energy
* `rama`:           ramachandran preferences
* `omega`:          omega dihedral in the backbone
* `total`:          total energy for structure/residue

Evaluating results
------------------
The score application outputs one summary file (`scorefile.sc`), containing summary information for all input files, and then one scored PDB file per input structure.
In the summary file, there is one line per input structure, listing the total energy from each term over the whole conformation.
At the end of each of the output PDB structures is a table with a per-residue breakdown of score terms.
The table is useful for identifying poorly scoring residues or analyzing the contribution of different score terms to the total for a particular residue or residues.
