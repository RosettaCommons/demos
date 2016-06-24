# RNA threading/mutation

KEYWORDS: NUCLEIC_ACIDS STRUCTURE_PREDICTION RNA

# Author
Rhiju Das, rhiju@stanford.edu

## Brief Description

Take an RNA template for homology modeling, and thread in a new sequence.

## Abstract

RNA homology modeling often involves taking a template coordinate file and threading on a new sequence. In the midst of the 2011 "RNA puzzles" community-wide blind trials, I hacked together some Rosetta code to do this.


# Running

```
$> $ROSETTA3/bin/rna_thread.default.linuxgccrelease -in:file:fasta rosetta_inputs/3bo3_REARRANGE_to_GIR1.fasta -s rosetta_inputs/3bo3_REARRANGE.pdb  -o 3bo3_REARRANGE_to_GIR1_thread.pdb -seq_offset 63
```

(where `$ROSETTA3`=path-to-Rosetta/main/source)

Note that insertions will *not* be modeled -- there will just be a chainbreak at those residues, and the residue numbering will skip.

If you don't have any alignment gaps in the target model, you can skip the creation of the alignment.fasta file, and just supply -seq <target sequence> on command line. See demo for rna_mutate.


