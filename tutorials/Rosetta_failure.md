## Limitation of Rosetta and things that can go wrong

Rosetta is a powerful tool that has been shown to be successful in many applications including many described either in tutorials and demos or the published paper. However, as with many other softwares, it has its own limitations. This part of the tutorial is meant to give you an idea of what some of these limitations are. Overall, one should bear in mind that no computational software is perfect and not all the things Rosetta outputs will work experimentally. In fact, user inputs have huge effects on iproving the outputs. Moreover, computational design should be used as a guide to reduce experimental efforts tather than a written-in-stone final solution. 

##### Scoring and scoring-dependent biases

The current Rosetta score, talaris 2014, has its own flaws specially when it comes to electrostatic interactions. //Asked Hahnboem to give me a paragraph on this.

In many of the examples in the tutorials, we sort our results based on the final score. You should note that in real applications, addition, you should note that when you only choose based on the final score, you may obtain structures that are not optimal. One classic example of this is the "all ALA" helix designs. Another example is the design of structures with many aromatic rings or with highly charged surfaces. This is because Rosetta tries to maximize the final score by adding this favorable interactions which are not possible in real world, or are not desirable for later experimental analysis (highly charged surfaces make crystallography hard, Ala residues cannot provide any special interactions despite being highly favored for helical regions).

*what the user can do*

The best thing to do is to be aware of the limitations of scoring function in Rosetta and to always, ALWAYS, inspect your structures carefully. Some of the biases mentioned can be controlled by newly added [aa composition](https://www.rosettacommons.org/docs/latest/rosetta_basics/scoring/AACompositionEnergy) functionlaity. Also, try to sort your structure not just by the final score but also by some additional filter values that you set. Take a look at other columns in the scorefile and see if they look normal. Check [score tutorial](scoring) for more information on the score function.

##### Rosetta's defaults

One thing that is not immediately obvious to new users is some of Rosetta defaults that are no longer in use but are still default. You can find some of them in demos and tutorials being set to false in the options files (no-OptH as an example). Another example is that when you make a new peptide, for example using PeptideStubMover, rosetta's default fot omega is 0, while the actual value should be 180. These hidden defaults can affect the outputs a lot without you knowing what was wrong. Another good example can be found in the [fold tree tutorial](fold_tree.md). You can see that how using a default tree can result in non-realistic outputs. One other example is that the constraint scores are OFF by default so if you have constraints, you need to turn them on. 

*what the user can do*

Rosetta developers are working constantly to update the settings so that it works better, however some of these issues are inevitable. Again, always check your output structures and see if it MAKES SENSE. Do the bonds and angles and torsions look right? Are the amino acids in the correct ramachandran bin? If you are not familiar with a task, read [the demos and tutorials](Demos, Tutorials, and Protocol Captures), check what options they set and read the [documentation](https://www.rosettacommons.org/docs/latest/Home) for more information on each. 

##### Fragments and Rosetta

One of the major ways Rosetta work is through fragment picking. While Rosetta score function has actual calculations, it also have some statistical components to it that biases the designs ot the native structures (see [above]()). So, if you are designing a topology that has never been observed in PDB or very rarely, you are in risk of not being able to score it well or fold it. Another subsequent result of using fragments for de nove structure predicition is that if you are trying to predict strcuture of something that has no close strcutures in the PDB, there is a good chance you WONT'T be able to get the correct structure.
Also, if you use a noncanonical residue in the file, scoring may be problematic (see [below](#Soemthing_old_something_new:noncanonicals).

*what the user can do*

For issues like this, you need to find other approaches to the problem.

##### The curse of oversampling

One thing to remember is that "you get what you ask for". At the end, Rosetta performs what you ask it to do. You don't want to bias your results by a series of assumptions and restrcitions. For example if you run a [MonteCarlo](link to Monte Carlo) algorithm that is set to optimize harshly on interface energy, you may end up getting designs with very good interface energy but with unreasonable number of hydrophobics, or unreasonable orientations. 

*what the user can do**

Try to avoid oversampling. Never bias your search by assumptions you have about the problem. Run the task with different settings. Check the outputs and look carefully.

##### Something old something new:noncanonicals

Rosetta was originally designed for proteins. Also, the folding prediction and scoring is highly dependent on PDB fragments. In recent years more and more new functionalities has been added to Rosetta. Now it can fairly well handle nucleic acids in the structures. With given params file and constraints, and sometimes slight changes in databse, ligands and metals can be handled. D-amino acids can be used relatively easily and many noncanonical amino acids and theri [rotemer libraries](Optimizing_Sidechains_The_Packer) are available either with the main package or upon request. Glycan handling and relaxation is on the way. However, there are still many challenges facing noncanonical usage in Rosetta, particularly in terms of scoring and analysis. Beta amino acids cannot be scored well, many noncanonical amino acids are lacking, and we still lack a good way of scoring structures longer than 40 aa that contain D amino acids.

*what the user can do*

Many researchers are working hard to increase Rosetta's power to be able to handle and score these noncanonicals. But if you are using a noncanonical, you need to familiarize yourself with the current limitations and latest advances. You may have to remove the noncanonical region from part of your scoring analysis or use other methods to rank your results.

##### Final points

**Rosetta cannot do magic**. It is just a software with many functionalities designed to be used by users to guid in solving their problems. YOU, as the users, are part of the process that cannot be ignored. In fact, in all the steps of the process the user has a significant part: guiding the software by controlling inputs and outputs (check [input and output tutorial](input_and_output)), defining what needs to be done and what parts of the structure we want to change through task operations and residue selectors (check (rosetta scripts tutorial)(rosetta_scripting)), and at the end analyzing the data. Always check the results not just the parts you are intereted in but the whole structure. See what is happening and always ask: **Does this make sense?**

And if you have any questions, you can always rely on the [rosetta community](https://www.rosettacommons.org/support).
