#Core Rosetta Concepts

KEYWORDS: CORE_CONCEPTS GENERAL

Author: Frank Teets(teetsf@gmail.com), Parisa Hosseinzadeh (parisah@uw.edu)

[[_TOC_]]

##Rosetta Numbering

Rosetta numbering of residues in a protein structure, sometimes called Pose numbering, is distinct from the residue numbering in the input PDB files ("PDB numbering"). Rosetta numbering always begins at 1 for the first residue and increases by one for each residue, ignoring chain designation. This has some implications for residue identification. Heuristically, when Rosetta reads in or outputs residue numbers, the presence of a chain identifier indicates PDB numbering, whereas the absence of a chain letter normally indicates Rosetta numbering. If Rosetta complains about residues that do not exist in the PDB, it may be using the other numbering scheme.

##Pose 

Rosetta's internal representation of a protein structure is called a Pose. It contains the atoms in the structure and the chemical connectivity as well as other information like how Rosetta interconverts internal/Cartesian coordinates (the FoldTree), the current energies, any constraints applied to the structure, and annotations from the original PDB. Specific protocols can add more information to the Pose as well.

##Internal Coordinates

Rosetta stores the structure in both Cartesian and internal coordinate representations, and handles the interconversion between the two. This allows protocols to do more efficient sampling of torsion space, manipulating internal coordinates directly. This has two implications that should generally be kept in mind:

1. Any change to a backbone torsion angle without a compensatory change to another torsion angle will result in a lever-arm effect, with residues down-tree of the move moving in 3D-space, even if they were marked as not moving.
2. By default, Rosetta does not optimize (or score) bond lengths and angles, either keeping the input geometries or building them from "ideal" values. Specialized protocols to optimize and score bonding geometry exist, but must be explicitly invoked to correct bad bond geometries.

##Monte Carlo Sampling

Rosetta is stochastic software; given the huge search space available, we walk through it randomly rather than exhaustively. Accordingly, a simulation must be run many times to produce a population of results, called *decoys*, which may be analyzed statistically. Generally, the number of output structures can be controlled by the `-nstruct` (number of structures) option. See [Rosetta at Different Scales](https://www.rosettacommons.org/docs/latest/getting_started/Rosetta-on-different-scales) for a general sense of the recommended number of output structures for different tasks, and refer to the documentation and demos for details on specific protocols.

##Rotamers

Rotameric isomers, or rotamers, indicate the position of a given amino acid side chain; they consist of an amino acid identity and the chi angles of the side chain heavy atoms. Rosetta uses the [Dunbrack](http://dunbrack.fccc.edu/) backbone dependent rotamer library, which is derived from high quality protein crystal structures.

Rosetta does not draw a fundamental distinction between sampling within and between amino acid rotamer populations; that is to say, both packing and design are handled by choosing side chains from a library, with the difference between packing and design determined by which amino acids are used to populate the library.

##Limitations of Rosetta and things that can go wrong

Rosetta is a powerful tool that has been shown to be successful in many applications, including many described either in tutorials, demos and published papers. However, as with most cutting-edge scientific software, it has its limitations. This part of the tutorial is meant to give you an idea of what some of Rosetta's limitations are and how to minimize their impact. One should always bear in mind the programming maxim "garbage in, garbage out": the quality of Rosetta's output is limited by the care in which input data is made, and how well your modeling task matches the assumptions made by the protocol used. Moreover, computational design and modeling should be viewed as a tool for reducing the search space of a problem, not providing a single, unequivocally "right" solution.  

##### Scoring and scoring-dependent biases

The current Rosetta score, ref2015, has certain known limitations (see the list below). While these are not major, they can occasionally result in certain behaviors that can affect the quality of the final output. 

- preference for aromatics (over aliphatic or sulfur-containing (=MET))
- preference for hydrogen bonding over other electrostatic interactions
- weak electrostatic repulsion term related to above mentioned terms 
- weak penalty on solvation of buried polars 
- poor description of context dependence of electrostatic interaction
- preference for helices (can be observed from structure prediction)
- poor description of pre-Pro / cis-Pro preferences
- wrong rotameric preferences in certain AAs (this is somewhat app dependent, but can be poor for MET, ASP, ASN,TRP)

In many of the examples in the tutorials, we sort our results based on the final score. These scores are general, and a naive sort by score alone may return structures that are not optimal due to the above mentioned limitations in the score function. One classic example of this is the "all ALA" helix designs. Another example is the design of structures with many aromatic rings or with highly charged surfaces. This is because Rosetta tries to maximize the final score by adding this favorable interactions which are not possible in real world, or are not desirable for later experimental analysis (highly charged surfaces make crystallography hard, Ala residues cannot provide any special interactions despite being highly favored for helical regions).

*What You Can Do:*

The best thing to do is to be aware of the limitations of scoring function in Rosetta and to always, ALWAYS, inspect your structures carefully. Some of the biases mentioned can be controlled by newly added [aa composition](https://www.rosettacommons.org/docs/latest/rosetta_basics/scoring/AACompositionEnergy) functionality. Also, try to sort your structure not just by the final score but also by some additional filter values that you set. Take a look at other columns in the scorefile and see if they look normal. See the [[scoring tutorial|scoring]] for more information on the score function.

##### Rosetta's defaults

One thing that is not immediately obvious to new users is some Rosetta defaults are no longer recommended but are still default. You can find some of them in demos and tutorials being set to false in the options files (no-OptH as an example). Another example is that when you make a new peptide, for example using PeptideStubMover, Rosetta's default for omega is 0, while the actual value should be 180. These hidden defaults can affect the outputs a lot without you knowing what was wrong. Another good example can be found in the [[fold tree tutorial|fold_tree]]. You can see how using a default tree can result in non-realistic outputs. One other example is that the constraint scores are OFF by default so if you have constraints, you need to turn them on. Some functions do this automatically, but not all. 

*What You Can Do:*

Rosetta developers are working constantly to update the settings so that it works better, however some of these issues are inevitable. Again, always check your output structures and see if it MAKES SENSE. Do the bonds and angles and torsions look right? Are the amino acids in the correct Ramachandran bin? If you are not familiar with a task, read [[the demos and tutorials|Home]], check what options they set and read the [documentation](https://www.rosettacommons.org/docs/latest/Home) for more information on each. 

##### Fragments and Rosetta

Rosetta relies heavily on data derived from experimental structures - both for statistical potentials and for sampling through fragment picking and rotamer libraries. This means that Rosetta results are biased toward structural motifs seen in existing experimental structures. So if you are designing or modeling a topology or structure that has never or only very rarely been observed in the PDB, you are at risk of not being able to model it accurately. This also means that chemical entities which are found only rarely in the PDB, like noncanonical amino acid, may not be sampled or scored properly (see [below](#limitations-of-rosetta-and-things-that-can-go-wrong_something-old-something-new-noncanonicals).

*What You Can Do:*

For issues like this, fragment-free algorithms are available.

##### The curse of oversampling

One thing to remember is that "you get what you ask for". At the end, Rosetta performs what you ask it to do. For example if you run a [Monte Carlo](#monte-carlo-sampling) algorithm that is set to optimize strictly on the interface energy, you may end up getting designs with very good interface energy but with unreasonable number of hydrophobic amino acids, or unreasonable sidechain orientations. 

*What You Can Do:*

Try to avoid oversampling and do not optimize on other score terms to the detriment of the actual score.

##### Something old something new:noncanonicals

Rosetta was originally designed for proteins, and the folding prediction and scoring is highly dependent on PDB fragments. In recent years more and more new functionalities has been added to Rosetta. Now it can fairly well handle nucleic acids in the structures. With given params file and constraints, and sometimes slight changes in database, ligands and metals can be handled. D-amino acids can be used relatively easily and many noncanonical amino acids and their [[rotamer libraries|Optimizing_Sidechains_The_Packer]] are available either with the main package or upon request. Glycan handling and relaxation is on the way. However, there are still many challenges facing noncanonical usage in Rosetta, particularly in terms of scoring and analysis. Beta amino acids cannot be scored well, many noncanonical amino acids are lacking, and we still lack a good way of scoring structures longer than 40 aa that contain D amino acids.

*What You Can Do*

Many researchers are working hard to increase Rosetta's power to be able to handle and score these noncanonicals. But if you are using a noncanonical, you need to familiarize yourself with the current limitations and latest advances. You may have to remove the noncanonical region from part of your scoring analysis or use other methods to rank your results.

##### Final points

Rosetta has a limited ability to detect nonsensical starting conditions and a series of internal sanity checks. It has no way of knowing the applicability of its results to a particular problem. You, the end user, must evaluate both your inputs and Rosetta's results on their scientific merits, both statistically and experimentally. 

And if you have any questions, you can always rely on the [Rosetta community](https://www.rosettacommons.org/support).
