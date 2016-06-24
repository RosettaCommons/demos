Protein DNA Interface Design
==============================

KEYWORDS: NUCLEIC_ACIDS INTERFACES DNA

This docuemnt is last updated in 2016 by Parisa Hosseinzadeh.

In this demo, we provided several simple scripts for handling protei\_DNa interfaces and esigning them. For more information, please refer to [rosettaDNA documentation](manual/applications/app_rosettaDNA.dox).

*Step 1*: Preparing Structures
Rosetta should be capable of reading in any protein-DNA complex in the modern PDB format
Some heavy-atom nucleotide residues may need to be converted by hand into canonical nucleotides. (e.g.: Iodocytosine...)
Alternative conformations will be ignored by Rosetta. In this demo, you will be using the pdb file 2h7h.pdb, provided in starting_files directory. Copy the file to your directory:
```
  $> cp starting_files/2h7h.pdb .
```

In order to obtain the initial score, run the command below:

(where `$ROSETTA3`=path-to-Rosetta/main/source)
```
  $> <$ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease @score.options > score.log
```
Now, pack the structure by running the script below. This packing routine finds optimal sidechain configurations for the wildtype amino acids in the vicinity of the target nucleotide base pair(s) (which are specified in the pack.script file as "dna_defs"):
```
  $> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease @pack.options > pack.log
```
Note the changes between the scores before and after packing. Compare the pack_2h7h and score_2h7h structures. If you do not have the files, we have provided them for you in the output directory. Please note that these are only simple short examples and for real designs, you need more sampling. You also need to change the residue numbers based on your specific application.

*Step 2*: Design Protein

This routine designs low energy amino acid identities and sidechains in the vicinity of the target nucleotides. Run the script below:
```
  $> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease @design.options > design.log
```
Example of output structures and scores are found in output directory.

*Step 3* Multi-state Design

This routine optimizes a starting population of different single-state designs for specificity toward the target DNA sequence vs. its competitors, using a genetic algorithm that involves mutation and recombination of the most sequence-specific solutions over multiple generations. Run this script as an example:
```
  $> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease @multistate.options > multistate.log
```
Example of output structures and scores are found in output directory.

Typically, this multi-state genetic algorithm routine should be run with a population >= 1000 protein sequences, for >= 50 generations.
