#Optimizing Side-Chains: Introduction To The Packer
Tutorial by Vikram K. Mulligan, Ph.D.
Created 20 June 2016.

[[_TOC_]]

## Summary

Rosetta's primary algorithm for optimizing side-chains is called the <i>packer</i>.  By the end of this tutorial, you should understand:
- The types of problems that the packer solves.
- The way in which the algorithm works, and what you can do to ensure optimal packer behaviour.
- How to invoke the packer using the <i>fixbb</i> application.
- How to control the packer's behaviour.
- What applications and algorithms invoke the packer.

The tutorial also introduces:
- Common ways of invoking the packer from the RosettaScripts scripting language.
- Common ways of controlling the packer from the RosettaScripts scripting language.

## The Problem of Optimizing Side-Chains

A common task in Rosetta is the optimization of side-chains.  We might think about the problem as follows: let's suppose that we have a structure (which we will call the <i>pose</i>) with a fixed backbone conformation.  At each position in the structure, we have a list of discrete possibilities (which we term <i>rotamers</i>) for the side-chain, and we would like to select one rotamer for each position such that the combination of rotamers represents the lowest-energy solution.  This is the problem solved by the packer.

This problem is actually quite a difficult one: given N possibilities at each position in an M-residue protein, there are N to the power of M possibilities.  This rapidly becomes an astronomical number of possibilities -- for example, 3 rotamers at each of 100 positions would be about 5x10<sup>47</sup> possible combinations.  This makes exhastive enumeration impossible.  To solve this problem, the packer uses Monte Carlo methods (discussed in detail further on).  This means that the packer is stochastic, and that the solution returned, while likely to be a <i>good</i> solution, is not guaranteed to be the <i>best</i> solution.

> <b>Repeated packer runs are likely to yield a variety of similar solutions near the global optimum; none of these will necessarily <i>be</i> the best possible solution.</b>

## Invoking the packer through the <i>fixbb</i> application

To illustrate what the packer does, let's run it.  The simplest way to run the packer is by running the Rosetta <i>fixbb</i> application: its sole <i>raison d'être</i> is to call the packer.  Navigate to the demos/public/fixbb directory, and run the following:

```
<path_to_Rosetta_directory>/main/source/bin/fixbb.default.linuxgccrelease -in:file:s 1l2y.pdb -in:file:fullatom -resfile resfile.txt -nstruct 5 >log.txt 2>err.txt &
```

You may need to change "linuxgccrelease", in the above, to whatever is appropriate given your operating system and compiler.

This application packs the side-chains of the input structure (the trp cage mini-protein, 1l2y.pdb).  Five output structures, from five separate runs, are produced.  If you compare these structures to the input structure, you'll find that Rosetta chooses slightly different rotamers for the side-chains, as compared to the input.  This is to be expected, particularly given the discrete nature of rotamers: the truly "best" rotamer might lie between two rotamers tested, and may never be sampled.

Work through the [[rest of the demo|../../public/fixbb/Readme.md]].  This teaches about how the packer can be tweaked, both at the commandline and with configuration files called <i>resfiles</i>, to control the amount of sampling, the time taken for a run, and the likelihood of converging to the optimal solution.

## The Sequence Design Problem

It is important to note that the packer problem, as described above, makes no assumptions about the nature of the candidate side-chains at each position.  The lists of possibilities can just as easily contain different side-chain identities as it can contain different conformations of the same side-chain.  The packer is therefore a powerful tool (and, indeed, the primary tool in Rosetta) for designing amino acid sequences.

To see the packer design a sequence, navigate to the demos/public/fixbb_design directory and run the following:

```
<path_to_Rosetta_directory>/main/source/bin/fixbb.default.linuxgccrelease -in:file:s 1l2y.pdb >log.txt &
```

This will produce the output files 1l2y_0001.pdb and score.sc.  If you open 1l2y_0001.pdb in a PDB viewer and compare it to the input file, 1l2y.pdb, you'll see that the sequence has changed considerably.  The rotamers chosen, however, should be interacting with one another reasonably favourably -- that is, there shouldn't be side-chains occupying the same space (clashing), for example.  Note that the only flag here is the one to specify our input file; that is, we're not passing any options to the fixbb application in this case to control the behaviour of the packer.  This brings up a very important point:

> <b>The default behaviour of the packer is to <i>design</i> at every position, allowing every rotamer of each of the 20 canonical amino acids.</b>

If you complete the rest of the fixbb_design demo, you'll learn how to use resfiles to control the behaviour of the packer, allowing only certain positions to be designed, and only with certain amino acid residue types.
