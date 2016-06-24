# Fixed backbone design of RNA

KEYWORDS: NUCLEIC_ACIDS DESIGN RNA 

This code is intended to carry out fixed backbone design of RNA sequences given an input backbone.

This demo redesigns a 'UUCG' tetraloop on a single-base pair RNA 'helix', as a small 6-nucleotide test case. As illustration, only 3 designs are output. It takes about 15 seconds to run. Run:(where `$ROSETTA3`=path-to-Rosetta/main/source)

```
$> $ROSETTA3/bin/rna_design.default.linuxgccrelease @flags
```

The output will show up in:

```
 chunk001_uucg_RNA.pack.txt 
```

with scores in

```
 chunk001_uucg_RNA.pack.out
```

and PDBs in S_0001.pdb, S_0002.pdb, etc.  
The typical sequence output is cuuggg (native is cuucgg). 
