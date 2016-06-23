# Stepwise Enumerative Assembly

KEYWORDS: STRUCTURE_PREDICTION GENERAL

## Author
Rhiju Das, rhiju@stanford.edu
## Protocol Name
Was called: "StepWiseAssembly", but that was too generic. 
We need a catchier name. SEA? Enumerosetta?

## Brief Description

An enumerative ansatz for high resolution RNA and Protein Folding
Rhiju Das, Parin Sripakdeevong, Das Lab, Stanford Biochemistry
Tuesday Aug. 3, 2010

## Abstract

High-resolution structure modeling is severely limited by difficulties in conformational sampling. Current Rosetta approaches use a low-resolution sampling stage, fragments from experimental structures, or a Monte-Carlo-like search -- typically all three. We describe an alternative "ansatz" that builds well-packed models in small steps, enumerating several million conformations for each residue, and covering all possible paths. We present results on non-canonical RNA motifs as well as highly irregular protein loops that have been intractable for prior fragment assembly or analytic loop closure approaches. In all cases, the method either reaches atomic accuracy or exposes flaws in Rosetta’s high-resolution energy function. Blind tests on a tetraloop-receptor motif are underway, as well as extension of the method to more complex systems such as aptamers and knotted cyclotides.


## Running
### Example Rosetta Command Line:

```
  stepwise_protein_test.macosgccrelease  -rebuild -fasta mini_1alc.fasta  -cluster:radius 0.1  -score:weights score12_no_hb_env_dep.wts  -pack_weights pack_no_hb_env_dep.wts -add_peptide_plane -align_pdb mini_1alc_H.pdb -native mini_1alc_H.pdb  -s1 noloop_mini_1alc_H.pdb  -input_res1 `seq 1 11` `seq 20 28` -sample_res 12 -out:file:silent build_first_residue.out  -calc_rms_res `seq 12 19` -fixed_res `seq 1 11` `seq 20 28` -database ~/minirosetta_database/
```

### Example Overall Command Line (if overall protocol is run via a script or other program)
```
  grinder_dagman.py  -loop_start_pdb noloop_mini_1alc_H.pdb  -loop_res  `seq 12 19`   -align_pdb mini_1alc_H.pdb  -fasta mini_1alc.fasta  -nstruct 200  -native mini_1alc_H.pdb   -final_number 50 -denovo 1

[This is available in input/README_SETUP. ]

[The "seq" phrases make use of linux's seq command; for macs you either need to install this, or type out 12 13 14 15 16 17 18 19 for `seq 12 19`.]
```
The result should be:
```
  protein_build.dag 
```
An example is available in output/.

"dag" stands for "directed acyclic graph". This is in the format recognized by CONDOR's "dagman", which was the original queuing paradigm.  However, I found it really slow; lots of ltency in queuing jobs. So I ended up writing my own scripts to run on Stanford's BioX2 cluster, which uses LSF (Load Sharing Facility):
```
  SWA_pseudo_dagman_continuous.py  -j 20 protein_build.dag  

[This is available in input/README_SUB]
```
We have also written scripts to carry out the calculation on condor and torque clusters, and are almost finished with versions that use Amazon's EC2/S3. It usually takes a day or so to rewrite what we have for arbitrary systems. Our next step may be to write a version for Amazon's ElasticMapReduce, but it gets complicated since we actually have a series of map/reduce steps, not just one. 


## Versions

1. The right Rosetta version:

https://svn.rosettacommons.org/source/branches/das_lab/mini
Revision: 36561

[This branched off trunk in the winter of 2009.]

2. Several scripts to generate a loop modeling job are in:
 https://svn.rosettacommons.org/source/workspaces/rhiju/python
 Revision: 36561

[the scripts  include:
  grind_dagman.py
 stepwise_post_process_cluster.py
 stepwise_post_process_combine_and_filter_outfiles.py
 stepwise_pre_process_setup_dirs.py
 extract_lowscore_decoys.py

 and various helper python scripts...
]

 Also note that some of the python scripts look in ~rhiju for other python scripts -- I think if you change the paths at the top of grind_dagman.py to rosetta and to the python directory, you'll be good to go. If someone fixes this to be more generic, please let me know. We have to do it before publication anyway...

3. Finally, the queuing scripts are in:

 https://svn.rosettacommons.org/source/branches/das_lab/SWA_dagman_python 
 Revision: 36561

  I should probably combine the important scripts from #2 with the scripts in #3.

## References to published works using this protocol

This is unpublished (we're waiting for blind predictions to be tested... almost there!). Let us know if you extend the work, so we can trade tips and queuing strategies.

## Other Comments:
1. The run can be sped up if you use a subset of the protein around the loop of interest... 

2. To start extending the method to design, I think you just have to comment out some lines in src/protocols/swa/protein/StepWiseProteinPacker.cc:
    ```
    234		pack_task_->restrict_to_repacking();
    313		pack_task_->restrict_to_repacking();
    ```
    but obviously there will be more steps to calibrate the procedure for design applications.

3. We also have extensive scripts written for RNA building, (stepwise_rna_test.cc), and are planning to refactor and unify RNA/protein into one code flow in the Winter of 2010. 
