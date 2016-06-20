#Optimizing Side-Chains: Introduction To The Packer
Tutorial by Vikram K. Mulligan, Ph.D.
Created 20 June 2016.

[[_TOC_]]

## Summary

Rosetta's primary algorithm for optimizing side-chains is called the <i>packer</i>.  By the end of this tutorial, you will understand:
- The types of problem that the packer solves.
- The way in which the algorithm works, and what you can do to ensure optimal packer behaviour.
- How to invoke the packer using the <i>fixed_bb</i> application.
- How to control the packer's behaviour.
- What applications and algorithms invoke the packer.

The tutorial also introduces:
- Common ways of invoking the packer from the RosettaScripts scripting language.
- Common ways of controlling the packer from the RosettaScripts scripting language.

## The Problem of Optimizing Side-Chains

A common task in Rosetta is the optimization of side-chains.  We might think about the problem as follows: let's suppose that we have a structure (which we will call the <i>pose</i>) with a fixed backbone conformation.  At each position in the structure, we have a list of discrete possibilities for the side-chain.
