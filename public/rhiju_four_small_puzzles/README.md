# Four Small Puzzles
KEYWORDS: CORE_CONCEPTS GENERAL
## Author
Rhiju Das, rhiju@stanford.edu

## Brief Description

  Follow up to talk:

  An enumerative ansatz for high resolution RNA and Protein Folding
  Rhiju Das, Parin Sripakdeevong, Das Lab, Stanford Biochemistry
  Tuesday Aug. 3, 2010
  
  associated with paper: "Four Small Puzzles that Rosetta Doesn't Solve",
  submitted in Jan. 2011 to RosettaCon 2010 Special Collection of 
  PLoS One

## Abstract

A complete macromolecule modeling package must be able to solve the simplest structure prediction problems. Despite recent successes in high resolution structure modeling and design, the Rosetta software suite fares poorly on deceptively small protein and RNA puzzles, some as small as four residues. To illustrate these problems, this manuscript presents extensive Rosetta results for four well-defined test cases: the 20-residue mini-protein Trp cage, an even smaller disulfide-stabilized conotoxin, the reactive loop of a serine protease inhibitor, and a UUCG RNA tetraloop. In contrast to previous Rosetta studies, several lines of evidence indicate that conformational sampling is not the major bottleneck in modeling these small systems. Instead, approximations and omissions in the Rosetta all-atom energy function currently preclude discriminating experimentally observed conformations from de novo models at atomic resolution. These molecular “puzzles” should serve as useful model systems for developers wishing to make foundational improvements to this powerful modeling suite.


## Running
# Example Rosetta Command Line

 A full set of READMEs is available in input_files. These are directly copied from my run directories, so they better work.

## Example Overall Command Line (if overall protocol is run via a script or other program)

 A few runs involve 'StepWise Assembly' and require a master python script to setup the job, and another master python script to queue up the resulting computation, which is described as a directed acyclic graph. A full set of README_SETUPs and README_SUBs is available in input_files. These are directly copied from my run directories, so they better work too.  The runs work on an LSF queuing system, as is available on Stanford's BioX2 cluster;  they should also run fine via Condor's DAGMAN, though that queuing system has a lot of latency. [We have also written scripts to carry out the calculation on condor and torque clusters, and are almost finished with versions that use Amazon's EC2/S3. It usually takes a day or so to rewrite what we have for arbitrary systems. Our next step may be to write a version for Amazon's ElasticMapReduce, but it gets complicated since we actually have a series of map/reduce steps, not just one. ]


## Versions

1. The right Rosetta version:
```
https://svn.rosettacommons.org/source/branches/das_lab/mini
https://svn.rosettacommons.org/source/branches/das_lab/minirosetta_database
Revision: 36561

[This branched off trunk in the winter of 2009.]
```
2. Several scripts to generate a loop modeling job are in:
```
 https://svn.rosettacommons.org/source/workspaces/rhiju/python
 Revision: 40199

[the scripts  include:
  grind_dagman.py
 stepwise_post_process_cluster.py
 stepwise_post_process_combine_and_filter_outfiles.py
 stepwise_pre_process_setup_dirs.py
 extract_lowscore_decoys.py

 and various helper python scripts...
]
```

3. Finally, the queuing scripts are in:
```
 https://svn.rosettacommons.org/source/branches/das_lab/SWA_dagman_python 
 Revision: 40199
```
## References 

The protein and RNA de novo modeling make use of published protocols, e.g.:

i. Das R, Qian B, Raman S, Vernon R, Thompson J, et al. (2007) Structure prediction for CASP7 targets using extensive all-atom refinement with Rosetta@home. Proteins 69: 118-128.

ii. Mandell DJ, Coutsias EA, Kortemme T (2009) Sub-angstrom accuracy in protein loop reconstruction by robotics-inspired conformational sampling. Nat Methods 6: 551-552.

iii. Das R, Karanicolas J, Baker D (2010) Atomic accuracy in predicting and designing noncanonical RNA structure. Nat Methods 7: 291-294.

There are two manuscripts in preparation on StepWise Assembly for proteins and RNA.

## Other Comments
