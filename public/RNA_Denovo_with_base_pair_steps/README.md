# Fragment assembly of RNA (FARNA) with base pair steps (-bps_moves).

README written in Nov. 2015 by Rhiju Das.

This demo recovers the sarcin-ricin loop 1Q9A, leaving out first base pair, which is opened up in xtal structure.

It uses -bps_moves, knowledge-based fragments that include two base pairs separated by at most a 3-nt bulge. 
(At the time of demo release, there are no benchmarks or publications for this functionality.)

Make sure you have rna_tools set up (https://www.rosettacommons.org/docs/latest/application_documentation/rna/rna-tools), since
this uses rna_denovo_setup.py.
```bash
export ROSETTA_TOOLS=<path/to/Rosetta/tools>
$> source $ROSETTA_TOOLS/rna_tools/INSTALL
```

Go into input_files/ 

```bash
$> cd input_files/
```

and run

```bash
source README_SETUP
source README_FARFAR
```

To run an extremely fast version (for testing purposes), run

```bash
$> source README_SETUP.short
$> source README_FARFAR
```
