Structure prediction of an RNA-protein complex
=====================================================================================

KEYWORDS: NUCLEIC_ACIDS RNA INTERFACES

Written in April 2018 by Kalli Kappel (kappel at stanford dot edu).  

This demo shows how to model the 3D structure of an RNA-protein complex starting from a protein structure and RNA sequence.


## Setting up the demo: 

1. Install Rosetta RNA tools. See instructions and documentation [here](https://www.rosettacommons.org/docs/latest/application_documentation/rna/RNA-tools).

## Brief explanation of input files: 

1. `fasta.txt` is a file specifying the sequence of the RNA-protein complex.
2. `secstruct.txt` is a text file specifying the secondary structure for the RNA-protein complex in dot-bracket notation. Secondary structure for the protein should be specified by dots. The secondary structure should be the same length as the sequence found in the fasta file. 
3. `unbound_protein.pdb` is the unbound protein structure (you could alternatively use a homology model of the protein or the bound protein structure if it's known). This will be treated as a rigid body throughout the run. This structure must contain all of the protein residues that are in your `fasta.txt` file (this protocol does not build protein residues from scratch).
4. `RNA_helix.pdb` is a structure of a 3-residue helix (generated with `rna_helix.py` in RNA tools). This will be treated as a rigid body throughout the run.
5. As an alternative to `unbound_protein.pdb` and `RNA_helix.pdb`, `unbound_protein_and_RNA.pdb` is a structure containing both the unbound protein structure and the RNA helix. Putting both structures in the same file will fix the relative rigid-body orientation of the two and no docking moves will be performed.
6. `flags` contains all the flags for the run (or alternatively `flags_no_dock` for the run where the RNA helix is fixed relative to the protein).

## Running the demo:

Models of the RNA-protein complex will be built with the Rosetta fold-and-dock method, which combines FARNA RNA folding with RNA-protein docking. First, type:

```
rna_denovo @flags
```

This will take several minutes to run and will generate 5 structures. For a normal run, it is typically best to generate several thousand structures.  

The final structures will be found in the silent file `2qux_fold_and_dock.out`.

To extract PDB files from the silent file, type:

```
extract_lowscore_decoys.py 2qux_fold_and_dock.out 5
```

This will create 5 PDB files named 2qux_fold_and_dock.out.1.pdb, 2qux_fold_and_dock.out.2.pdb, etc.

Alternatively, to build models of the RNA-protein complex keeping the rigid-body orientation of the RNA helix and the protein fixed, first type:

```
rna_denovo @flags_no_dock
```

Again, this will take several minutes to run and will generate 5 structures in a silent file named `2qux_fold_and_dock_fix_rigid.out`. To extract PDB files from the silent file, type:

```
extract_lowscore_decoys.py 2qux_fold_and_dock_fix_rigid.out 5
```

These structures should be viewed in a molecular graphics program of your choice, e.g. Pymol.

For reference, example output for this demo is provided in the `example_output` directory. 

## Additional information

See additional documentation [here](https://www.rosettacommons.org/docs/latest/application_documentation/rna/rnp-modeling).
