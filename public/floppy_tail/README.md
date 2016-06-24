Floppy Tail Demo
================

KEYWORDS: STRUCTURE_PREDICTION GENERAL

This demo contains the starting structure for the original published use of 
FloppyTail, and is a copy of FloppyTail's integration test.  To run the demo 
(in the inputs subfolder):

```bash
$> cd inputs  
$> $ROSETTA3/FloppyTail.linuxclangrelease @options
```

Read the options file first to set local paths as needed.

Please refer to FloppyTail's documentation (in the [online manual](https://www.rosettacommons.org/docs/latest/floppy-tail.html), or in your copy of the code at 
Rosetta/documentation/application_documentation/structure_prediction/floppy-tail.md, and its [publication](http://www.ncbi.nlm.nih.gov/pubmed/19945379) 
Kleiger G, Saha A, Lewis S, Kuhlman B, Deshaies RJ. Rapid E2-E3 assembly and 
disassembly enable processive ubiquitylation of cullin-RING ubiquitin ligase 
substrates. Cell. 2009 Nov 25;139(5):957-68. PubMed PMID: 19945379.
Note that the paper is primarily about biology, not modeling.

Note that the provided outputs are from the integration test, which runs in ~30 
s.  You will need to run FloppyTail for much longer (both longer individual 
runs, and many trajectories) to get scientifically useful results.  The options 
file and documentation provide details on how to do that.
