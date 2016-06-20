Full Atom Representation vs Centroid Representation
===================================================

Need for a reduced represenation
--------------------------------
In an ideal world with infinite time and computer power, we could perform all of our simulations with all atoms. In practice, trying to perform extensive backbone sampling while including all side chain atoms is impractical at best.

The first problem is that having all atoms is expensive, because calculating interactions between all atom pairs grows rapidly (n<sup>2</sup>) with the number of atoms. The bigger problem is that fully atomic conformational space is _very rugged_, so most moves are rejected by Monte Carlo.

Centroid representation
-----------------------

To get around this problem, poses are often converted into centroid mode for portions of a protocol that require extensive sampling (for example, the initial stages of [ab initio structure prediction](https://www.rosettacommons.org/demos/latest/public/abinitio/README)). In centroid mode, the backbone remains fully atomic, but the representation of each side chain is simplified to a single _pseudo-atom_ of varying size. For protein backbones, this representation preserves five backbone atoms for each amino acid: nitrogen (N), the alpha carbon (CA), the carbonyl carbon (C), the carbonyl oxygen (O), and the polar hydrogen on nitrogen. The side chain is replaced by the `CEN` atom whose radius and properties (polarity, charge, etc.) are determined by the residue's identity.

Centroid score functions
------------------------
Centroid score functions are kind of vague, in the same way that the protein representation is kind of fuzzy. This has a disadvantage in terms of interpreting their results, but a huge advantage in that the energy landscape is not nearly as rugged, and sampling very different conformations is easier.

Converting to full atom
-----------------------
After large-scale sampling in centroid mode, poses are generally converted back to their all-atom representation for refinement, which generally entails some combination of side chain repacking and minimization. This allows Rosetta to more accurately score interactions between side chains and other finer details of the protein's structure.

Demo
----
[[images/1qys.png]]
[[images/1qys_centroid.png]]