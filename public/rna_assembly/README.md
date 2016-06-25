# RNA Assembly

KEYWORDS: STRUCTURE_PREDICTION NUCLEIC_ACIDS RNA

# Author

This README was written in Sep. 2011, by Rhiju Das (rhiju@stanford.edu); updated in Feb. 2012 after directory restructuring. Thanks to P. Kerpedjiev for suggestions.

Updated in Jun. 2016 by Kristin Blacklock (kristin.blacklock@rutgers.edu).

# Brief Info
This demo illustrates a protocol to assemble models of large RNAs by first building their helical stems and inter-helical motifs, and then putting them together.

It is being published in a (primarily experimental) paper "A two-dimensional mutate-and-map strategy for non-coding RNA structure" by W. Kladwang, C. VanLang, P. Cordero, and R. Das (2011), Nature Chemistry.

# Running the demo

The example input files are in rosetta_input; you may wish to copy them locally with the command:

```
$> cp rosetta_inputs/* .
```

Before running the command to set up the rna assembly job, make sure you have the correct paths to the rna_tool/bin directory in your environment variables. This can be done by running,
```bash
   export ROSETTA_TOOLS=<path/to/Rosetta/tools>
   source $ROSETTA_TOOLS/rna_tools/INSTALL
$> export RNA_TOOLS=$ROSETTA_TOOLS/rna_tools
$> python $RNA_TOOLS/sym_link.py
$> export PATH=$PATH:$RNA_TOOLS/bin/
$> export PYTHONPATH=$PYTHONPATH:$RNA_TOOLS/bin/
```

Everything needed to run the job is created by the command:

```
    $> ./scripts/setup_rna_assembly_jobs.py  add.fasta add_secstruct.txt 1y26_RNA.pdb add_mutate_map_threetertiarycontacts.cst exe_extension=".default.linuxgccrelease"
```

The first two arguments are required -- the sequence_file and the secondary structure file [either in dot/bracket notation, or specifying Watson/Crick base pairs as pairs of numbers]. 

The last three arguments are optional; they supply the native pdb and any constraints, here derived from a high-throughput "mutate-and-map" strategy for RNA structure determination. The exe_extension= option specifies the executable extension to use for the Rosetta executables.

You may need to change the path to your rosetta executable ('EXE_DIR') in setup_rna_assembly_jobs.py [in which case you should get a warning!]. The script currently assumes that you are in the rosetta_demos/public/rna_assembly directory within the rosetta codebase.



Then run the Rosetta commands in :

```
    README_STEMS
    README_MOTIFS
    README_ASSEMBLE
```

To run the extremely short version of this tutorial (mostly for testing purposes), run
```bash
$> ./README_STEMS
$> ./README_MOTIFS.short
$> ./README_ASSEMBLE.short
```

You can see examples of these files and their output in example_output/. Please note that for these files I changed the 'nstruct' commands to create 100 models per motif. In reality you will want to make 2000-4000 MOTIF models, and then several thousand ASSEMBLE models. [You can just use one STEM model per helix, as that is supposed to be an ideal helix.] For some scripts to generate lots of models on a computer cluster, see note below.

The final 'outfile' is  add_assemble.out. If you don't have the file ready, you can copy:
```
$> cp ./example_output/add_assemble.out .
```

We can extract models from it using:(where `$ROSETTA3`=path-to-Rosetta/main/source)

```
$> $ROSETTA3/bin/extract_pdbs.default.linuxgccrelease -in:file:silent add_assemble.out -tags S_000001
```

or using scripts like my `extract_lowscore_decoys.py`.

Caveat: The above protocol is a bit inflexible in that the motifs are modeled separately from each other -- if a loop/loop interaction occurs in the final global model it will not really be modeled correctly by the isolated loops. We are working on iterative methods to tackled this global de novo assembly question. For now the protocol seems to work well if there are experimental constraints that e.g., connect the loops.

# Appendix: Generating lots of models on a cluster

We often use either condor or LSF ('bsub') queuing systems on clusters, and have some handy scripts to set up these jobs, which are included in the Rosetta distribution.

The syntax (e.g., for README_MOTIFS) is:

```
    python <path/to/Rosetta/tools>/rna_tools/bin/rosetta_submit.py  README_MOTIFS MOTIF_OUTPUT 50
```

this will create the directory MOTIF_OUTPUT with 50 subdirectories 0/, 1/, etc. for 50 parallel jobs. You can then run:

```
    source bsubMINI
```

or

```
    condor_submit condorMINI
```

to kick off jobs. The rosetta_submit.py script is pretty straightforward and should be easy to edit for any kind of cluster queuing system.

After the jobs have been running for a while, or have finished, you can type:

```
    python <path/to/Rosetta/tools>/rna_tools/bin/easy_cat.py MOTIF_OUTPUT 
```

to concatenate your files.  [Do a similar thing with README_ASSEMBLE after getting the results of MOTIF].
