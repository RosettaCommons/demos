Anchored Design
==============

KEYWORDS: INTERFACES DESIGN

 This document describes how to use the AnchoredDesign protocol, both in benchmarking and design mode.  As the protocol's components [AnchorFinder](http://www.rosettacommons.org/docs/latest/anchor-finder.html), [AnchoredPDBCreator](http://www.rosettacommons.org/docs/latest/anchored-pdb-creator.html), and [AnchoredDesign](http://www.rosettacommons.org/docs/latest/anchored-design.html) are reasonably extensively documented elsewhere, this protocol capture is meant to be used alongside that online documentation. Presented at RosettaCon2010 (in poster form) was a description of the protocol itself, plus benchmarking results, plus some early design results.  The accompanying paper ([Lewis SM, Kuhlman BA. Anchored design of protein-protein interfaces. PLoS One. 2011;6(6):e20872. Epub 2011 Jun 17.](http://www.ncbi.nlm.nih.gov/pubmed/21698112) (pubmed link)) describes only benchmarking results, but the tools to do design are described here.  A paper on design results is forthcoming.

Note that this protocol capture is somewhat focused on just the AnchorFinder portion (the least important part of the process), because the other portions are documented elsewhere but AnchorFinder largely is not.  


Contained here:
* Instructions on choosing appropriate benchmarks – AnchorFinder or otherwise
* Instructions on preparing those structures – for benchmarking or via AnchoredPDBCreator
* How to benchmark AnchoredDesign against those structures
* How to use AnchoredDesign to design interfaces
* command lines/option files, and discussion of the options

Not contained here:
* A speck of code – that lives in your Rosetta release.
* Submission scripts for running jobs.  I don't know your architecture.  It's all MPI compatible so it's not hard.

Sort-of contained here:
* The raw_documentation directory includes copies of the doxygen-style manual documentation for AnchorFinder, AnchoredDesign and AnchoredPDBCreator.  These copies are guaranteed NOT to be up to date; look in your copy of the code or the web links above instead.

Overview
--------
The purpose of the AnchoredDesign protocol is to create new protein-protein interactions using information (the anchor) borrowed from a known interaction at the same interface of one partner.  Because this protocol is intended to design protein-protein interactions, the obvious test is to see whether it can recover known structures of such interactions.  This document is massively overconcerned with benchmarking because it accompanies the paper in which AnchoredDesign is benchmarked; unless you are actually trying to replicate the benchmarking you can ignore most of those details and skip to the design tools.

The protocol modifies a loop region around an anchor in designing binders.   Selection of structures for benchmarking therefore requires interface loops with anchor regions.

The ideal anchor has several qualities:
* one or a few contiguous residues – the protocol can only have one anchor
* many interactions across an interface – this is the whole point of an anchor
* not embedded in the middle of a secondary structure element (helix/strand) – the anchor needs to be in a loop because the interface space will be sampled via loop remodeling of the anchor loop

Examples of anchors might be a phosphotyrosine inserting into an SH2 domain, a polyproline sequence binding an SH3 domain, etc.

For design, one would choose an anchor based on one's target.  For benchmarking, you are free to choose anything that has an interface loop with a good anchor.  To select benchmarking structures, I wrote the AnchorFinder protocol.  AnchorFinder has some value in highlighting which residues might make good anchors for a given target (although computational alanine scanning, not covered here, is more likely to be useful).

After finding suitable structures with the help of AnchorFinder, the next step is to pick anchors and loops out of those structures, and in general prepare them for Rosetta's use (removing solvent atoms, etc).  At this point you're ready to run AnchoredDesign.


Compiling AnchorFinder
----------------------
(This section applies to benchmarking only)

Your Rosetta code distribution should include an application called AnchorFinder. If you wish to search large numbers of PDBs for potential anchors (I searched a local copy of the entire PDB structure set), then you will wish to modify the code slightly before running it.  Running any part of Rosetta against huge numbers of unprepared, straight-from-the-PDB structures is challenging because the PDB reader in Rosetta is not robust against nonstandard file formats, etc.

To compile AnchorFinder such that it will be robust, examine the manual documentation on RobustRosetta (also included in this protocol capture).  Briefly, this documentation describes changes that A) make Rosetta slower (thus they aren't on by default) and B) cause it to throw C++ exceptions when it hits errors instead of crashing.  The job distributor catches the errors, skips the bad structures, and continues.  You must recompile after making these changes.

You do not want to use compiled executeables OTHER than AnchorFinder with these changes made – they will significantly slow the code down.  AnchorFinder is quite fast so it's not a problem.

When running AnchorFinder, watch your memory usage.  When I used it, there was a patch in the JobDistributor which deleted starting poses for PDBs that had already been processed.  This patch was rejected by the community and since been replaced by a different patch to do the same thing; AnchorFinder is a run-once sort of thing so it has not been tested against the new method.

Using AnchorFinder
------------------
(This section is minimally relevant if not benchmarking)

At this point you should have a compiled copy of AnchorFinder with the necessary changes to the code.  You can then list your PDBs in one or many -l files (or -s) for use in Rosetta.  The format for Rosetta's -l flag is one path per line:

    A.pdb
    B.pdb
    C.pdb
    ...

Depending on your available architecture, it may be better to split the run up into many -l on separate processors.  I don't know what's best for you.

If you want to do the whole PDB – it's a good idea to skip the largest PDB files ahead of time, particularly ribosome structures.  These take a very long time to process through the PDB reader, and due to heavy nucleic acid content are skipped anyway.  You can also either toss the NMR structures ahead of time or use the -obey_ENDMDL flag to only read the first model.

AnchorFinder will automatically remove nonprotein atoms from the Poses before examination.  It also skips anything that is monomeric, has no protein residues, or smaller than 20 residues after processing.

It will then look through the structures searching for regions with certain command-line-defined characteristics.  These characters are:
* length of windows for consideration - 4 or 5 or 6, etc, contiguous residues at a time?  This flag is -window_size.  I suggest 5 residue windows.
* What fraction of this window should have loop secondary structure as assigned by DSSP?  -loopness controls this.  It takes it as a decimal between 0-1, I suggest 0.6 (which translates to 3/5 residues for a 5 residue window) to 1 (all residues loop).
* How many cross-interface interactions per residue should the window have?  I suggest a minimum of 4.  This translates to 20 (redundancy included) cross-interface interactions for a 5 residue window.  By redundancy, I mean residues 43 and 44 on chain A can both interact with residue 234 on chain B and it will count as two interactions.  Specify this with -nbrs_per_residue.
* What file name should the good interactions be printed to?  I leave it as an exercise to the reader to pick their own file name.  Specified with -bestoutfile; defaults to goodfile.out.

Running AnchorFinder, while not particularly slow, is still something you only want to do once.  The defaults suggested above produce lots of output, which can then be further processed quickly without reloading PDBs.  To expedite this, AnchorFinder produces two levels of output.  All residues have their data printed to a file named (pdbname).data – you can reprocess this to get data for differing window lengths, loopnesses, etc.  Windows passing the loopness and interactions filters are printed to the specified output file.

A suggested options file for AnchorFinder is available with this document.

Interpreting AnchorFinder Results
---------------------------------
(This section is minimally relevant if not benchmarking)

After you've run AnchorFinder, you'll have a fairly large pile of output: pdbname.data for all pdbs, plus goodfile.out for the better windows.

`pdbname.data` looks like this:

Rows are residues, columns are chains, data are neighbors in that chain for each residue

    residue chain   PDBdata DSSP    1       2
    1       1       2 D     L       7       0
    2       1       3 D     L       10      0
    3       1       4 D     L       14      0
    ...

The columns are residue and chain in Rosetta numbering, residue/chain in PDB numbering, DSSP value, and then N columns for the N chains in the protein.  The number in those columns is the number of cross-interface neighbors on that chain for that position.

`goodfile.out` looks like this:

    PDB pdb2vk1 window 45 loopness 5 nbrs 0 28 0 0 start 46 A pymol select pdb2vk1 and chain A and resi 46-50
    PDB pdb2vk1 window 108 loopness 5 nbrs 0 25 0 0 start 109 A pymol select pdb2vk1 and chain A and resi 109-113
    PDB pdb2vk1 window 109 loopness 5 nbrs 0 36 0 0 start 110 A pymol select pdb2vk1 and chain A and resi 110-114
    PDB pdb2vk1 window 110 loopness 5 nbrs 0 46 0 0 start 111 A pymol select pdb2vk1 and chain A and resi 111-115
    PDB pdb2vk1 window 111 loopness 5 nbrs 0 46 0 0 start 112 A pymol select pdb2vk1 and chain A and resi 112-116
    PDB pdb2vk1 window 112 loopness 5 nbrs 0 47 0 0 start 113 A pymol select pdb2vk1 and chain A and resi 113-117

Each line identifies the PDB, the window number, its loopness, its number of neighbors on each chain in the PDB (variable # of columns), the starting residue PDB numbering for the window, and a Pymol selection for the window.

Inputs and outputs for this stage from a convenience sample (PDBs 3cy?) are included with this protocol capture.

At this point, the data is yours to play with.  I searched for windows with large numbers of neighbors on only one chain using sifter.py (included), then sorted for those with the largest number of neighbors (sort -n -k1 `input`).  After that it was all manual filtering to choose structures for the benchmarks.

Choosing Loop and Anchor — Benchmarking
---------------------------------------
(This section applies to benchmarking only)

OK, so you ran AnchorFinder, looked at the results, and/or picked what protein you want to run through AnchoredDesign.  How do you choose a loop/anchor?

If you ran AnchorFinder, look at the AnchorFinder result lines that came up as good:

    92 PDB pdb1zr0 window 526 loopness 5 nbrs 0 0 92 0 start 13 D pymol select pdb1zr0 and chain D and resi 13-17

Load this PDB into pymol (1zr0.pdb) and activate the suggested selection.  You'll see that it is in a surface loop of one partner which sticks an arginine straight into its binding partner – a perfect anchor.  (This is a chosen example; not all AnchorFinder hits are this nice.)

Choosing the anchor is entirely up to human effort; here the arginine 15 is an obvious choice.

For choosing loops, I just traveled up and down the chain in both directions until I hit secondary structure, significant backbone-backbone hbonding, or the protein core.  Here I'd choose a loop of D10 to L17 – more N-terminal than that affects the core, and more C-terminal affects a sheet.

Anchor and loop file specifications are included in the release documentation and the examples here.

Note that for the included example, the PDB has been renumbered from 1.  Scripts to do this are occasionally included with Rosetta distributions and not included here.  It will be convenient to also remove waters, ligands, etc.

If you are doing benchmarking, skip to the [Running AnchorDesign](#Running 
AnchorDesign) section.

Choosing a System and Anchor — Design
-------------------------------------
(This section applies to the design case only)

In the design case, you will be choosing your proteins based on what you want designed.  Your target is forced by what targets:
* have crystal structures, and
* are related to your biological problem.

Choosing an anchor then requires:
* a cocrystal of your target with some partner from which to source the anchor.

You can run this cocrystal through AnchorFinder and let it suggest anchors to you, but for one structure you can just look at it yourself.  Look for loops on the partner that insert into the target, or do computational alanine scanning, or examine the literature for mutations that disrupt the interface.

Choosing a Scaffold, Design Positions, and Loops — Design
---------------------------------------------------------
(This section applies to the design case only)

In the design case, you will be replacing your target's partner with some new scaffold to form a mostly de novo interface.  Your scaffold must meet a few requirements:
* flexible, mutateable surface loops (for AnchoredDesign to modify)
* experimentally tractable (hey, your funeral if it's not)
* whatever other functionality you need for your desired design

The protocol was written with the fibronectin monobody scaffold in mind.

Choosing which loops are flexible is dependent on biological knowledge of the scaffold.  In fibronectin's case, many papers have been published establishing the mutability of the BC and FG loops.

Choosing which positions are designable is similarly dependent on your scaffold.  AnchoredDesign carries the assumption that the non-anchor loop positions are designable, and non-loop positions are not, but nothing in the code enforces that.  Use a resfile (documented with the manual) to specify which positions are designable.  The code will automatically prevent design of the anchor (you can turn that off).  The code will automatically prevent design of positions that are not close to either the interface or a flexible loop (you cannot turn that off), so take care in specifying designable positions on opposite faces of your protein.  Proximity is redetermined at each design opportunity so positions peripheral to the interface may not be designed regularly.

Choosing Loop Lengths and Anchor Position
-----------------------------------------
(This section applies to the design case only)

OK, so you know which scaffold to use, and which anchor, and which target.  You are ready to create your starting structure for AnchoredDesign, in which the anchor will be inserted into the scaffold, and the anchor will be aligned properly to the target, dragging the scaffold and target together.  The protocol used for this is called AnchoredPDBCreator; further details are below.

One important part of conformational space that AnchoredDesign cannot search is the space of loop lengths and anchor positions.  You may want to try, for a loop of length N, all combinations of loops of length N-3 to N+3, or even more for long loops.  As you are designing the loop to form an interface, there is no reason to believe its native length is particularly relevant.  You will have to do this searching at this stage: create starting structures for all loop lengths, run them all through AnchoredDesign, and pick off the best ones later.

Loops can be shortened directly by just deleting residues mid-loop before handing the scaffold to AnchoredPDBCreator – it can insert a 3 residue anchor into a 6 residue window, and close the gap.  Loop lengthening must be done externally.  One way to lengthen loops is to manually modify a PDB to contain enough residues in the loop (copy-and-paste a residue, renumber as necessary), then use the loop_modeling executeable's build_initial mode to close the loop.  Further instructions are included in their own folder in this packet.

A paired space is anchor placement space.  Besides choosing which anchor to use (try several), exactly where it is placed within a loop can vary.  For a loop of length 7, and an anchor of length 2, (assuming a flexible residue on each side), you have the following 4 choices:

    X = scaffold
    - = loop
    A = anchor
    X1234567X
    X-AA----X
    X--AA---X
    X---AA--X
    X----AA-X

Again, this space is not searched by AnchoredDesign and must be searched by trying all the inputs.

Using AnchorPdbCreator
----------------------
(This section applies to the design case only)

AnchoredPDBCreator is the protocol which assembles an anchor, scaffold, and target into a starting structure for AnchoredDesign.  Its code documentation is included in this packet.

Briefly, AnchoredPDBCreator takes as input 4 files:
* The target structure, as a PDB, with the partner removed
* The anchor structure, drawn from the cocrystal with the target, containing ONLY the residues being used as an anchor, as a PDB
* The scaffold structure, as a PDB, with loop residues added/deleted as desired
* A scaffold_loop specification, which declares which residues in the scaffold are flexible and where the anchor insertion should occur.

It is ABSOLUTELY VITAL to recognize that AnchoredPDBCreator does NOT produce interfaces, it only produces starting structures for AnchoredDesign.  It is entirely plausible that its structures will have the target and scaffold totally eclipsed.  This is fine, AnchoredDesign will fix it.

AnchoredPDBCreator's results should be interpreted by analyzing ONLY the closure of the anchored loop.  Use the result with the best loop geometry.  Loop geometry can be measured by examining the LoopAnalyzerMover output tagged to the end of result PDBs:

LoopAnalyzerMover: unweighted bonded terms and angles (in degrees)

    position phi_angle psi_angle omega_angle peptide_bond_C-N_distance rama_score omega_score dunbrack_score peptide_bond_score chainbreak_score
     pos phi_ang psi_ang omega_ang pbnd_dst    rama  omega_sc dbrack pbnd_sc   cbreak
      17  -106.8   175.8     178.2    1.322   0.998    0.0342   7.01   -2.68   0.0182
      18  -82.33   64.67    -178.5    1.329   0.211    0.0217   3.11   -3.42   0.0203
      19  -83.63   149.4     177.2    1.329   -1.07    0.0795      0   -3.43    0.584
      20  -75.25   171.1    -178.7    1.329  -0.264    0.0161  0.348   -3.43   0.0151
      21  -58.53  -42.95     174.6    1.329   -0.58     0.294      0   -3.43      2.7
      22  -76.02   159.9    -179.8    1.326  -0.811  0.000404   0.97   -3.45   0.0424
      23  -72.63   130.1     179.4    1.325   -1.29   0.00372   0.24   -3.46   0.0281
      24  -94.91   116.5     179.8    1.323   -1.21   0.00028  0.721   -3.45   0.0694
      25  -65.42   150.7     179.4    1.335   -1.58     0.004      0   -3.32     1.38
      26  -64.68   147.9     179.1    1.323   -1.45    0.0079   1.61   -3.32    0.211
      27  -56.44  -66.68      -180    1.329    1.34  8.08e-30   7.87   -3.43 2.37e-05
      28  -124.4  -56.48     177.6    1.329    2.08    0.0568  0.608   -3.43   0.0533
      29  -124.1   28.78    -177.7    1.264   0.341    0.0542   2.39    2.65     2.07
      30   81.57  -134.3    -176.4    1.329      20     0.126   5.06    2.65    0.128
      31  -112.9   147.2     172.7    1.318  -0.744     0.538  0.534   -3.35     1.38
    total_rama 15.9674
    total_omega 1.23676
    total_peptide_bond -38.3223
    total_chainbreak 8.70689
    total rama+omega+peptide bond+chainbreak -12.4113

    LAM_total -12.4113

In this particular example, position 29 is clearly problematic: the peptide bond distance is too short, as reported by the pbnd_dst, pbnd_sc, and cbreak columns.

You should be running AnchoredPDBCreator for at least 100 trajectories before choosing a starting structure.

Running AnchorDesign
--------------------
If you are benchmarking, the crystal structure of the complex is the appropriate input for AnchoredDesign.  If you are designing, the best result from AnchoredPDBCreator is your starting structure.

The input files for AnchoredDesign provide an example with 1zr0 for running AnchoredDesign.  It is a heterodimer so you can pretend it was AnchoredPDBCreator sourced if you want.  (You can also look in the AnchoredDesign integration test at test/integration/tests/AnchoredDesign for such an input).

* The anchor file specifies which residues form anchor.
* The PDB file is the pdb.
* loopsfile_extended is for benchmarking – the extension column is true, which tells AnchoredDesign to forget the starting loop conformation before sampling
* loopsfile_native tells AnchoredDesign NOT to forget the starting loop conformation – this is probably what you would use for design, although there is no reason you can't reject the starting loop conformation for design.  (This is cheating for benchmarking, but useful for making relaxed natives)
* options is the command line options file.  The active options are for a simple “does it run” test; parameters for a longer running “real” test are included.
* A resfile is necessary if you wish to design, and is only active in that option set.  

Interpreting AnchorDesign — Benchmarking
----------------------------------------
(This section applies to the benchmarking case only)

If you are duplicating the benchmarking results, you passed the rmsd flag. AnchoredDesign will have output a lot of RMSD values allowing you to determine the performance of the protocol against the structures you chose to benchmark.  The paper describes the score versus RMSD metrics used to determine quality (including the I_sup_bb_RMSD, ch2_CA_RMSD, and loop_CA_sup_RMSD.  The structures themselves don't really matter; you are ensuring that the low-scoring structures have low RMSD.

Interpreting AnchorDesign — Design
----------------------------------
(This section applies to the design case only)

In the design case, the other fields of the AnchoredDesign output come in to play. There are three classes of output:
* scorefunction terms
* LoopAnalyzerMover output,
* InterfaceAnaylzerMover output.

Generally, you should rank your structures according to total_score (the Rosetta scorefunction).  This tells you what Rosetta thinks is best.

Next, you use the LoopAnalyzerMover output (described above) and InterfaceAnalyzerMover output to determine which structures have flaws not caught by total_score.  Toss structures that those filters think have problems.  Pick the ones you think are best, order the DNA, and pray.  When it works great, feel free to send me kudos, citations, or money!

InterfaceAnalyzerMover
----------------------
InterfaceAnalyzerMover output looks like this:

    Residues missing H-bonds:
    Residue 	 Chain 	 Atom 
    38 	 A 	 NE2
    101 	 A 	 OE1
    248 	 A 	 O
    250 	 A 	 O
    344 	 B 	 N
    384 	 B 	 O
    477 	 B 	 O

    pymol-style selection for unstat hbond res 
    select start_5411_unsat, /start_5411//A/38+101+248+250+ + /start_5411//B/344+384+477+

    pymol-style selection for interface res 
    select start_5411_interface, /start_5411//A/31+32+33+34+35+36+37+38+39+40+41+54+56+57+59+60+61+62+64+65+66+92+95+98+99+100+101+102+103+106+194+195+224+225+226+227+228+229+230+247+248+249+250+251+252+253+265+ + /start_5411//B/314+315+316+317+318+319+320+321+322+323+324+337+339+340+342+343+344+345+347+348+349+350+375+378+379+381+382+383+384+385+386+389+473+476+477+478+479+480+481+482+488+506+507+508+509+510+511+512+513+

The first section documents where Rosetta thinks there are unsatisfied hydrogen bonds at the interface.  This code is known to be oversensitive to missing bonds, but it's better than nothing.

The next sections print PyMOL selections for interface residues for easier visualization.

InterfaceAnalyzerMover also includes columns into the scorefile:

    dSASA_int 2396.33
    dG_separated -35.3379
    dG_separated/dSASAx100 -1.47467
    delta_unsatHbonds 7
    packstat 0
    dG_cross -27.6963
    dG_cross/dSASAx100 -1.15578
    AllGly_dG -2.83564
    cen_dG -10.3844
    nres_int 96
    per_residue_energy_int -1.15006
    side1_score -361.478
    side2_score -267.353
    nres_all 520
    side1_normalized -1.25513
    side2_normalized -1.15238
    complex_normalized -1.74616
    hbond_E_fraction 0.368537

Most of these are experimental and not useful (and not part of AnchoredDesign; InterfaceAnalyzerMover has other clients).  The useful ones are dG_separated/dSASAx100, which measures the Rosetta energy of binding per unit area of SASA (scaled by a factor of 100).  This ensures you pick an interface that is energetic for its size, not large but sloppy.
