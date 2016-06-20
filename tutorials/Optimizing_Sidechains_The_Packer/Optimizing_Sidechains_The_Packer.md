#Optimizing Side-Chains: Introduction To The Packer
Tutorial by Vikram K. Mulligan, Ph.D.
Created 20 June 2016.

[[_TOC_]]

## Summary

Rosetta's primary algorithm for optimizing side-chains is called the <i>packer</i>.  By the end of this tutorial, you will understand:
- The types of problem that the packer solves.
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

To illustrate what the packer does, let's run it.  Navigate to the demos/fixbb directory, and run the following:

```
<path_to_Rosetta_directory>/main/source/bin/fixbb.default.linuxgccrelease -in:file:s 1l2y.pdb -in:file:fullatom -resfile resfile.txt -nstruct 5 >log.txt 2>err.txt &
```

## The Sequence Design Problem

It is important to note that the packer problem, as described above, makes no assumptions about the nature of the candidate side-chains at each position.  The lists of possibilities can just as easily contain different side-chain identities as it can contain different conformations of the same side-chain.  The packer is therefore a powerful tool (and, indeed, the primary tool in Rosetta) for designing amino acid sequences.
> <b>The default behaviour of the packer is to <i>design</i> at every position, allowing every rotamer of each of the 20 canonical amino acids.</b>

