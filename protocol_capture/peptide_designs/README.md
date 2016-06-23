Peptide Backbone and Sequence Design
====================================
KEYWORDS: DESIGN GENERAL
Author: Deanne Sammond

RosettaCon Talk:
* Computational design of a new protein-protein interface between Gi1 and a 
  redesigned RGS14 GoLoco, Deanne Sammond, Dustin Bosch, Glenn Butterfoss, 
  Mischa Machius, David Siderovski, Brian Kuhlman, Kuhlman lab, 2010, Session 3 
  on Wednesday August 4th

---

Our project is the computational design of a new high-affinity protein-protein 
interface.  Our model system is an x-ray crystal structure of  Gi1 bound to the 
GoLoco domain from the RGS14 protein.  RGS14 GoLoco spans two domains of Gi1, 
with the C-terminal random coil region binding to the all-helical domain of 
Gi1.  We removed this C-terminal portion of GoLoco, replacing the random coil 
with a de novo designed alpha helix.  The redesigned GoLoco binds to Gi1 with a 
dissociation constant of 810nM, the correct binding of the newly designed 
GoLoco was confirmed using disruptive mutations at the Gi1:GoLoco interface, 
and the correctness of the computational design was assessed with by x-ray 
crystallography.

This protocol builds (or extends) a backbone for a peptide bound to a target 
protein, then designs a low-energy sequence. 

Running the protocol
--------------------

#### Important flags:

    -ex1, -ex2, -exOH, -extrachi_cutoff 1 all seem to be very important for the sequence design run.

#### Example Rosetta Command Line:

    rosetta.mactel aa input_pdb _ -s g000.pdb -loops
    rosetta.mactel -design -l list_of_pdbs -tail -begin 342 -end 351 -chain_ -series bb -protein g000 -resfile g000_resfile -ex1 -ex2 -extrachi_cutoff 1 -exOH -no_his_his_pairE -tight_hb -try_both_his_tautomers 

### How to generate tail designs:

This protocol uses 2 separate rosetta runs — one is centroid mode to build 
backbone coordinates and the other is a design run to find a low-energy 
sequence — and 2 additional scripts.  Step-by-step instructions are below:

1. Generate fragments - you can do this using the [[Robetta 
   server|http://robetta.bakerlab.org/fragmentsubmit.jsp]]. I named my fragment 
   files `aag00003_04.200_v1_3` and `aag00009_04.200_v1_3`.

2. Make starting structure using createTemplate.pl and the g000.zones file.  I 
   also included a fasta file (see g000\_.fasta) because I felt that setting the 
   sequence improved the quality of the centroid models more effectively than 
   using constraints.  For example:

        createTemplate.pl -zonesfile g000.zones -fastafile g000_.fasta -parentpdb gpep_nat.pdb -outpdb g000.pdb

   The resulting file should look something like g000.pdb.  The side-chains are 
   removed from all sequence positions, and the region that will be redesigned 
   is removed.  NOTE: The input file has the sequence positions re-numbered in 
   "Rosetta numbers", so the original pdb starts with sequence position 30, 
   contains a "TER" between chain A and B, and chain B starts with sequence 
   position 496.  But my modified pdb (gpep_nat.pdb) starts with position 1, 
   with no "TER" between chains A and B, and the numbering is sequential 
   through B.

   input: g000.zones, g000_.fasta, gpep_nat.pdb  
   output: g000.pdb

3. Making centroid models

   * Copy fragments to working directory

   * Make a loop file to specity what residues can move.  See g000.loops as an 
     example.  In the loop file, the 1st numer is the # of positions to be 
     built, second number is sequence position to start design, with the final 
     number being the sequence position where the design will end. 

   * Make constraint file to make the designed region fall into the desired 
     location OR direct the redesigned peptide toward the desired orientation 
     using the starting sequence.  I did the latter.  In other words, part 2) 
     above can generate a pdb file with all alanines in the region that will 
     be designed OR the sequence can be specified with a fasta file (see 
     above) so that big hydrophobics fall in buried regions and hydrophilics 
     fall in solvent exposed regions, etc.  An example of a constraint file is 
     g000_.cst.

   * Copy over starting structure g000.pdb from step 2.

   * Run command line:

            rosetta.mactel aa input_pdb _ -s g000.pdb -loops 

   input: g000.loops, g000.cst (I didn't use constraints)  
   output example: aag000_0001.pdb

4. Merge centroid designs with fullatom starting structure (in this case 
   gpep1.pdb) using merge_pdb.csh like this:

        merge_pdb.csh gpep1_nat.pdb [list of pdbfiles].

   You will need to edit merge_pdb.csh if you want to change which residues are 
   being merged.  For example, if you start the design at sequence position 
   342, like we do here, check your gpep_nat.pdb file (original all-atom pdb 
   file) for the line # for the last atom in sequence position 341 and put this 
   # in after "head -", then check your centroid files for the first atom at 
   position 342 and put this # after "tail -".  This step is so that you don't 
   have to repack all of the gpep1.pdb positions during the fullatom 
   simulations.  (NOTE: When building with centroid mode, we don't use a TER in 
   between chains A and B.  The TER needs to be added back in.  Another merge 
   file can be used for this.)

   input: gpep_nat.pdb, list_of_pdbs (example of pdbs in list - aag000_0001.pdb)  
   output example: aag000_0001.m.pdb~ and with the TER added, aag000_0001.m.pdb

5. Making Fullatom models:

        rosetta.mactel -design -l list_of_pdbs (the *.m.pdb merged files from step 4 above WITH a TER added between chains A and B) -tail -begin 342 -end 351 -chain_ -series bb -protein g000 -resfile tail.resfile -ex1 -ex2 -extrachi_cutoff 1 -exOH -no_his_his_pairE -tight_hb -try_both_his_tautomers -linmem_ig 10 -output_hbond_info -decoystats -group_uns 

   * `-linmem_ig 10` is optional.  I used it because I was running on a 
     BlueGene and each node had very limited memory.  `-output_hbond_info`, 
     `-decoystats` and `-group_uns` are also optional.  I used those so the 
     output pdbs could be used with Ron Jacak's h-bond pymol plugin.

   * Interesting "feature" that seems to have appeared in this version - this 
     call to Rosetta is looking for fragment files named aag000l03_05.200_v1_3 
     and aag000l09_05.200_v1_3, whereas the fragment files I used when building 
     centroid models (3) were named aag00003_04.200_v1_3 and 
     aag00009_04.200_v1_3.  I just renamed the fragment files so I could get 
     this done quickly.

   input: g000_resfile  
   output example: aag000_0001.m_0001.pdb

See example of resfile (g000_resfile).  In this file the sequence positions of 
the designed region are allowed to vary, and any neighboring sequence positions 
on the target protein are allowed to relax.

Rosetta Version
---------------
SVN Revision: 29304  
https://svn.rosettacommons.org/source/trunk/rosetta++

