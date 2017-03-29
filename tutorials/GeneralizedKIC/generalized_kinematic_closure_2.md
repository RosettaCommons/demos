# Generalized Kinematic Closure Tutorial 2:
# Perturbing loops from a starting conformation with GeneralizedKIC
======================================

KEYWORDS: LOOPS SCRIPTING_INTERFACES

Tutorial by Vikram K. Mulligan (vmullig@uw.edu).  Created on 28 March 2017 for the Baker lab Rosetta Tutorial Series.

[[_TOC_]]

## Goals

At the end of this tutorial, you will understand:

- How to use the GeneralizedKIC mover to sample small perturbations of an existing loop.
- How to use the GeneralizedKIC mover with the GenericMonteCarlo mover to perform a Monte Carlo search of loop conformations.

## Perturbing loops with GeneralizedKIC

In the [[first tutorial|generalized_kinematic_closure_1.md]], we saw how one can construct an entirely new loop using the [PeptideStubMover](https://www.rosettacommons.org/docs/latest/PeptideStubMover), the [DeclareBond mover](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/DeclareBond), and GeneralizedKIC.  The GeneralizedKIC mover was used to sample entirely new conformations of the loop.  However, there are many situations in which one may wish to perturb a loop slightly from a starting conformation without randomizing its conformation completely.  This is particularly true when using GeneralizedKIC to explore variant conformations of a starting structure, or when one has used fragment insertion to build a loop crudely and wishes to refine the starting model.  In this tutorial, we will learn how to perturb loops with GeneralizedKIC.

## Exercise 2:  A RosettaScripts-Scripted Monte Carlo Search of Loop Conformational Space

### Inputs

For this exercise, we will be using one of the imperfect loop conformations from the first tutorial.  This simulates the design case, in which one may have built an initial, imperfect loop using a fragment-based or other method, and one wishes to improve one's initial model.

**The starting model for this tutorial:**
![The starting model for this tutorial](images/Exercise2_startingstruct.png)


## Conclusion

**TODO**

## Further Reading

Bhardwaj G, Mulligan VK, Bahl CD, Gilmore JM, Harvey PJ, Cheneval O, Buchko GW, Pulavarti SV, Kaas Q, Eletsky A, Huang PS, Johnsen WA, Greisen PJ, Rocklin GJ, Song Y, Linsky TW, Watkins A, Rettie SA, Xu X, Carter LP, Bonneau R, Olson JM, Coutsias E, Correnti CE, Szyperski T, Craik DJ, Baker D.  (2016).  Accurate de novo design of hyperstable constrained peptides.  _Nature_ 538(7625):329-335.

Mandell DJ, Coutsias EA, Kortemme T. (2009).  Sub-angstrom accuracy in protein loop reconstruction by robotics-inspired conformational sampling.  _Nat. Methods_ 6(8):551-2.

Coutsias EA, Seok C, Jacobson MP, Dill KA.  (2004).  A kinematic view of loop closure.  _J. Comput. Chem._ 25(4):510-28.

[GeneralizedKIC documentation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/composite_protocols/generalized_kic/GeneralizedKIC)

