# Relax a Large Structure

KEYWORDS: STRUCTURE_PREDICTION GENERAL

The FastRelax protocol minimizes the input structure according to Rosetta's force field.
The attached XML file (demo.xml) is an example rosetta script that specifies how to relax a part of a large protein.
The script defines one mover (FastRelax) with a MoveMap, which specified the residues to relax.
The mover is referenced in the PROTOCOLS section as the single operation.

An empty MoveMap selects all residues for relaxation.
The attached script forces the positions 1-800 and 1000-1200 to be fixed and allows the rest of the protein to be moved.
On 3E0C.pdb for example (a protein that has 1011 positions), this will only relax positions between 800-1000 (watch out: the numbering is according to the pose and might not correspond to pdb numbering if there are missing numbers in the pdb).
The bb and chi parameters specify if the backbone and/or side chains should be relaxed.
The repeats parameter control the length of the relax simulation.

Multiple residue spans can be specified and will be parsed in the specified order.
For example, the snippet below turns off relax for all residues, and then specified that only residue 100-200 should be relaxed:

```
	<Span begin=1 end=1011 chi=0 bb=0/>
	<Span begin=100 end=200 chi=1 bb=1/>
```

If not explicitly specified, the following parameter values are used for FastRelax:

```
	- scorefxn: score12
	- repeats: 8
	- task_operations: InitializeFromCommandline, IncludeCurrent, and RestrictToRepacking
```

Command line arguments to run this script: (where `$ROSETTA3`=path-to-Rosetta/main/source)

```
	$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease -s starting_files/3E0C.pdb -parser:protocol demo.xml -ignore_unrecognized_res
```

The flag ignore_unrecognized_res asks Rosetta to ignore any residue types it doesn't recognize (HOH in this example).


Additional flags that can help:

```
	-linmem_ig <int>
	might reduce the memory usage of packing (althought it didn't seem to help on this example). 

	-relax::min_type lbfgs_armijo_nonmonotone
	This uses lbfgs minimizer instead of the bfgs minimizer that uses less memory. You can control this amount of this saving by specifying (default=64): -optimization::lbfgs_M. If the protein has 200 positions, lbfgs with lbfgs_M set to 200 will behave exactly like bfgs. You can trade-off memory usage for running time by trying to minimize for a longer time. Use (default=200) -optimization::lbfgs_max_cycles
```

