# Fragment assembly of RNA (FARNA) with base pair steps (-bps_moves).
KEYWORDS: NUCLEIC_ACIDS RNA DENOVO STRUCTURE_PREDICTION
README written in Nov. 2015 by Rhiju Das.
Updated in July 2016 to only use `rna_denovo` app

This demo recovers the sarcin-ricin loop 1Q9A, leaving out first base pair, which is opened up in crytal structure.

It uses -bps_moves, knowledge-based fragments that include two base pairs separated by at most a 3-nt bulge. 
(At the time of demo release, there are no benchmarks or publications for this functionality.)

Go into input_files/ 

```bash
$> cd input_files/
```

and run

```bash
$> $ROSETTA3/bin/rna_denovo.default.linuxgccrelease @flags
```

