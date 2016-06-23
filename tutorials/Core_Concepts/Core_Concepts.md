#Core Rosetta Concepts
KEYWORDS: CORE_CONCEPTS GENERAL
Authors: Frank Teets(teetsf@gmail.com) and Parisah
Jun 23 2016
##Rosetta Numbering
Rosetta numbering of residues in a protein structure, sometimes called Pose numbering,  is distinct from the PDB numbering of the same; Rosetta numbering always begins at 1 and increases confluently according to the Fold Tree; functionally it behaves as though it were a PDB with a single chain beginning at residue 1. This has some implications for residue identification; inexplicable references to residues that do not exist in the PDB suggest that somewhere a residue is identified in the wrong mode.
##Pose 
Rosetta's internal representation of a protein structure is called a Pose. It contains the atoms in the structure, the atomic connectivity tree (the AtomTree), the energies scored from that structure, and the contents of the original PDB(if any). Specific protocols can add more information to the Pose as well.
##Internal Coordinates
Rosetta stores both cartesian coordinates and internal coordinates; with some exceptions, most protocols sample rotational space (and therefore manipulating internal coordinates directly). This has two implications that should generally be kept in mind:
1. any change to a backbone torsion angle without a compensatory change to another torsion angle will result in a lever-arm effect, with residues down-tree of the move moving in 3space, sometimes dramatically.
2. By default, Rosetta assumes that all covalent bond geometry is ideal; cartesian score terms are available to score bond lengths and angles, and should be used if scoring a structure that was designed or otherwise manipulated to potentially contain nonideal bond geometry.
##Monte Carlo Sampling
Rosetta is stochastic software; given the huge search space available, we walk through it randomly rather than exhaustively. Accordingly, a simulation must be run many times to produce a population of results, called *decoys*, which may be analyzed statistically. The number of results recommended for a given application is referred to by the name of the option controlling it, *nstruct*; see [Rosetta at Different Scales] for a general sense of these recommendations and refer to the documentation and demos for specifics.
##Rotamers
Rotameric isomers, or rotamers, indicate the position of a given amino acid side chain; they consist of an amino acid identity and the chi angles of the side chain heavy atoms. Rosetta does not draw a fundamental distinction between sampling within and between amino acid rotamer populations; that is to say, it both packs and designs side chains by choosing rotamers from a library, with the difference between packing and design determined by which amino acids are used to populate the library.

## Limitations of Rosetta and things that can go wrong

Rosetta is a powerful tool that has been shown to be successful in many applications including many described either in tutorials and demos or the published paper. However, as with any other software, it has its neither omnipotent nore infallible. This part of the tutorial is meant to give you an idea of what some of Rosetta's limitations are and how to minimize their impact. One should always bear in mind the programming maxim "garbage in, garbage out"; that is to say, Rosetta's output can only ever be as relevant as the provided assumptions. Moreover, computational design is a tool for reducing the search space of a problem, not providing a single unequivocally "right" solution.  

##### Scoring and scoring-dependent biases

The current Rosetta score, talaris 2014, has its own flaws particularly when it comes to electrostatic interactions. //Asked Hahnboem to give me a paragraph on this.

In many of the examples in the tutorials, we sort our results based on the final score. These scores are general, and a naive sort by score alone may return structures that are not optimal. One classic example of this is the "all ALA" helix designs. Another example is the design of structures with many aromatic rings or with highly charged surfaces. This is because Rosetta tries to maximize the final score by adding this favorable interactions which are not possible in real world, or are not desirable for later experimental analysis (highly charged surfaces make crystallography hard, Ala residues cannot provide any special interactions despite being highly favored for helical regions).

*What You Can Do:*

The best thing to do is to be aware of the limitations of scoring function in Rosetta and to always, ALWAYS, inspect your structures carefully. Some of the biases mentioned can be controlled by newly added [aa composition](https://www.rosettacommons.org/docs/latest/rosetta_basics/scoring/AACompositionEnergy) functionality. Also, try to sort your structure not just by the final score but also by some additional filter values that you set. Take a look at other columns in the scorefile and see if they look normal. Check [score tutorial](scoring) for more information on the score function.

##### Rosetta's defaults

One thing that is not immediately obvious to new users is some of Rosetta defaults that are no longer in use but are still default. You can find some of them in demos and tutorials being set to false in the options files (no-OptH as an example). Another example is that when you make a new peptide, for example using PeptideStubMover, rosetta's default fot omega is 0, while the actual value should be 180. These hidden defaults can affect the outputs a lot without you knowing what was wrong. Another good example can be found in the [fold tree tutorial](fold_tree.md). You can see that how using a default tree can result in non-realistic outputs. One other example is that the constraint scores are OFF by default so if you have constraints, you need to turn them on. Some functions do this automatically, but not all. 

*What You Can Do:*

Rosetta developers are working constantly to update the settings so that it works better, however some of these issues are inevitable. Again, always check your output structures and see if it MAKES SENSE. Do the bonds and angles and torsions look right? Are the amino acids in the correct ramachandran bin? If you are not familiar with a task, read [the demos and tutorials](Demos, Tutorials, and Protocol Captures), check what options they set and read the [documentation](https://www.rosettacommons.org/docs/latest/Home) for more information on each. 

##### Fragments and Rosetta

One of the major ways Rosetta work is through fragment picking. While Rosetta score function has actual calculations, it also have some statistical components to it that biases the designs ot the native structures (see [above]()). So, if you are designing a topology that has never been observed in PDB or very rarely, you are in risk of not being able to score it well or fold it. Another subsequent result of using fragments for de nove structure predicition is that if you are trying to predict strcuture of something that has no close strcutures in the PDB, there is a good chance you WONT'T be able to get the correct structure.
Also, if you use a noncanonical residue in the file, scoring may be problematic (see [below](#Soemthing_old_something_new:noncanonicals).

*What You Can Do:*

For issues like this, fragment-free algorithms are available.

##### The curse of oversampling

One thing to remember is that "you get what you ask for". At the end, Rosetta performs what you ask it to do. For example if you run a [MonteCarlo](link to Monte Carlo) algorithm that is set to optimize harshly on interface energy, you may end up getting designs with very good interface energy but with unreasonable number of hydrophobics, or unreasonable orientations. 

*What You Can Do:*

Try to avoid oversampling and do not optimize on other score terms to the detriment of the actual score.

##### Something old something new:noncanonicals

Rosetta was originally designed for proteins, and the folding prediction and scoring is highly dependent on PDB fragments. In recent years more and more new functionalities has been added to Rosetta. Now it can fairly well handle nucleic acids in the structures. With given params file and constraints, and sometimes slight changes in databse, ligands and metals can be handled. D-amino acids can be used relatively easily and many noncanonical amino acids and theri [rotemer libraries](Optimizing_Sidechains_The_Packer) are available either with the main package or upon request. Glycan handling and relaxation is on the way. However, there are still many challenges facing noncanonical usage in Rosetta, particularly in terms of scoring and analysis. Beta amino acids cannot be scored well, many noncanonical amino acids are lacking, and we still lack a good way of scoring structures longer than 40 aa that contain D amino acids.

*What You Can Do*

Many researchers are working hard to increase Rosetta's power to be able to handle and score these noncanonicals. But if you are using a noncanonical, you need to familiarize yourself with the current limitations and latest advances. You may have to remove the noncanonical region from part of your scoring analysis or use other methods to rank your results.

##### Final points

Rosetta has a limited ability to detect nonsensical starting conditions and a series of internal sanity checks. It has no way of knowing the applicability of its results to a particular problem. You, the end user, must evaluate both your inputs and Rosetta's results on their scientific merits, both statistically and experimentally. 

And if you have any questions, you can always rely on the [rosetta community](https://www.rosettacommons.org/support).
