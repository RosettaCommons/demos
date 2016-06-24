As one of the test to evaluate the design we dock the TS structure into the protein.  The docking is described in detail in Davis & Baker JMB 2009.  Here a new LG.params is needed as the zinc ion is kept in the protein and not free to move.  A constraint file is used for the zinc-ligand distance.  The number of alternative conformations to be sampled is controlled by -nstruct, which should normally be set considerably higher (thousands of samples are typical). It is set to 10 to allow the demo to complete quickly.
KEYWORDS: UTILITIES GENERAL
To run the docking application, use the following commandline:

```bash
<path_to_Rosetta_directory>/main/source/bin/ligand_dock.default.linuxgccrelease @flags
```

You may need to change ".default.linuxgccrelease" in the above to match your build, operating system, and compiler.  For example, if you were using the Intel icc compiler and a static debug-mode build, it would become ".static.linuxiccdebug".
