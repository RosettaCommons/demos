#Optimizing Side-Chains: Introduction To The Packer
Tutorial by Vikram K. Mulligan, Ph.D.
Created 20 June 2016.

[[_TOC_]]

## Summary

Rosetta's primary algorithm for optimizing side-chains is called the *packer*.  By the end of this tutorial, you should understand:
- The types of problems that the packer solves.
- The way in which the algorithm works, and what you can do to ensure optimal packer behaviour.
- How to invoke the packer using the *fixbb* application.
- How to control the packer's behaviour.
- What applications and algorithms invoke the packer.

The tutorial also introduces:
- Common ways of invoking the packer from the RosettaScripts scripting language.
- Common ways of controlling the packer from the RosettaScripts scripting language.

## The problem of optimizing side-chains

A common task in Rosetta is the optimization of side-chains.  We might think about the problem as follows: let's suppose that we have a structure (which we will call the *pose*) with a fixed backbone conformation.  At each position in the structure, we have a list of discrete possibilities (which we term *rotamers*) for the side-chain, and we would like to select one rotamer for each position such that the combination of rotamers represents the lowest-energy solution.  This is the problem solved by the packer.

This problem is actually quite a difficult one: given N possibilities at each position in an M-residue protein, there are N to the power of M possibilities.  This rapidly becomes an astronomical number of possibilities -- for example, 3 rotamers at each of 100 positions would be about 5x10<sup>47</sup> possible combinations.  This makes exhastive enumeration impossible.  To solve this problem, the packer uses Monte Carlo methods (discussed in detail further on).  This means that the packer is stochastic, and that the solution returned, while likely to be a *good* solution, is not guaranteed to be the *best* solution.

> **Repeated packer runs are likely to yield a variety of similar solutions near the global optimum; none of these will necessarily *be* the best possible solution.**

## Invoking the packer through the *fixbb* application

To illustrate what the packer does, let's run it.  The simplest way to run the packer is by running the Rosetta *fixbb* application: its sole *raison d'être* is to call the packer.  Navigate to the demos/public/fixbb directory, and run the following:

```
<path_to_Rosetta_directory>/main/source/bin/fixbb.default.linuxgccrelease -in:file:s 1l2y.pdb -in:file:fullatom -resfile resfile.txt -nstruct 5 >log.txt 2>err.txt &
```

You may need to change "linuxgccrelease", in the above, to whatever is appropriate given your operating system and compiler.

This application packs the side-chains of the input structure (the trp cage mini-protein, 1l2y.pdb).  Five output structures, from five separate runs, are produced.  If you compare these structures to the input structure, you'll find that Rosetta chooses slightly different rotamers for the side-chains, as compared to the input.  This is to be expected, particularly given the discrete nature of rotamers: the truly "best" rotamer might lie between two rotamers tested, and may never be sampled.

Work through the [[rest of the demo|../../public/fixbb/Readme.md]].  This teaches about how the packer can be tweaked, both at the commandline and with configuration files called *resfiles*, to control the amount of sampling, the time taken for a run, and the likelihood of converging to the optimal solution.

## The sequence design problem

It is important to note that the packer problem, as described above, makes no assumptions about the nature of the candidate side-chains at each position.  The lists of possibilities can just as easily contain different side-chain identities as it can contain different conformations of the same side-chain.  The packer is therefore a powerful tool (and, indeed, the primary tool in Rosetta) for designing amino acid sequences.

To see the packer design a sequence, open a terminal window, navigate to the demos/public/fixbb_design directory and run the following:

```
<path_to_Rosetta_directory>/main/source/bin/fixbb.default.linuxgccrelease -in:file:s 1l2y.pdb >log.txt &
```

This will produce the output files 1l2y_0001.pdb and score.sc.  If you open 1l2y_0001.pdb in a PDB viewer and compare it to the input file, 1l2y.pdb, you'll see that the sequence has changed considerably.  The rotamers chosen, however, should be interacting with one another reasonably favourably -- that is, there shouldn't be side-chains occupying the same space (clashing), for example.  Note that the only flag here is the one to specify our input file; that is, we're not passing any options to the fixbb application in this case to control the behaviour of the packer.  This brings up a very important point:

> **The default behaviour of the packer is to *design* at every position, allowing every rotamer of each of the 20 canonical amino acids.**

If you complete the rest of the fixbb_design demo, you'll learn how to use resfiles to control the behaviour of the packer, allowing only certain positions to be designed, and only with certain amino acid residue types.

## Working efficiently with the packer: TaskOperations

When we ran the fixbb application to repack side-chains without design, we saw that there was a tradeoff between speed, probability of convergence, and accuracy.  Including more rotamers often allows Rosetta to find better arrangements of side-chains, but at the cost of longer packer runs that might be less likely to reach their lowest-energy state.  Designing with the packer greatly increases the number of rotamers: instead of including rotamers for just one type of side-chain, the packer must consider one or more rotamers for *each* type of side-chain considered.  When packing, it can be a good idea to run short test-runs, paying attention to the number of rotamers being considered in the tracer output.  The first packing example, in which amino acid identities were fixed and only a small number of conformations were being considered at each position, produced the following in the ouput log, for example:

```
core.pack.pack_rotamers: built 256 rotamers at 20 positions.
```

With so few rotamers and such a low number of positions, the packing job runs almost instantly.  Later, when we turned on extra rotamers (still keeping the sequence fixed) the output message changed to:

```
core.pack.pack_rotamers: built 6021 rotamers at 20 positions.
```

This is still well within the range of doable for Rosetta, but at this point, packing will take some time (probably on the order of seconds per run).

Design, without enabling extra rotamers (which is what we did initially in the second demo), is also somewhat computationally expensive, but still manageable:

```
core.pack.pack_rotamers: built 4737 rotamers at 20 positions
```

However, we greatly simplify the problem and speed things up by restricting the packer to at most only a few choices at each position with a resfile, as we did in the second part of the second demo:

```
core.pack.pack_rotamers: built 410 rotamers at 16 positions.
```

In general, it is a good idea to limit the packer's options as much as possible in order to keep the search quick and convergent.  This brings up another important concept:

> **Each run of the packer can by controlled with one or more TaskOperations.**

TaskOperations can be passed to the packer in one of several ways, and different TaskOperations modify packer behaviour in different ways.  We've already seen one example: a user can control packer behaviour with a resfile specified with the "-resfile" flag on the commandline, which implicitly invokes the ReadResfile TaskOperation.  Resfiles can control amino acid identity, presence of extra rotamers, and other packer behaviours on a residue-by-residue basis, specified by residue index.  For more information about resfiles and their full features 
