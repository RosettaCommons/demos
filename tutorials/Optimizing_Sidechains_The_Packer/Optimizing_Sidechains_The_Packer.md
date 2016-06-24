#Optimizing Side-Chains: Introduction to the *Packer*

KEYWORDS: CORE_CONCEPTS GENERAL STRUCTURE_PREDICTION DESIGN   

Tutorial by Vikram K. Mulligan, Ph.D.   
Created 20 June 2016.

[[_TOC_]]

## Summary

Rosetta's primary algorithm for optimizing side-chains is called the *packer*.  By the end of this tutorial, you should understand:

- The types of problems that the packer solves.
- How to invoke the packer using the *fixbb* application.
- How to control the packer's behaviour.
- How to ensure optimal packer behaviour.
- What applications and algorithms invoke the packer.

The tutorial also introduces:

- Common ways of invoking the packer from the RosettaScripts scripting language.
- Common ways of controlling the packer from the RosettaScripts scripting language.
- A brief overview of how the packing algorithm works.

## The problem of optimizing side-chains

A common task in Rosetta is the optimization of side-chains.  We might think about the problem as follows: let's suppose that we have a structure (which we will call the *pose*) with a fixed backbone conformation.  At each position in the structure, we have a list of discrete possibilities for the side-chain, which we call *rotamers*, where a rotamer is a particular conformation of a particular residue type's side-chain.  We would like to select one rotamer for each position such that the combination of rotamers represents the lowest-energy solution.  This is the problem solved by the packer.

This problem is actually quite a difficult one: given N possibilities at each position in an M-residue protein, there are N to the power of M possibilities.  This rapidly becomes an astronomical number of possibilities -- for example, 3 rotamers at each of 100 positions would be about 5x10<sup>47</sup> possible combinations.  This makes exhaustive enumeration impossible.  To solve this problem, the packer uses Monte Carlo methods (discussed in detail further on).  This means that the packer is stochastic, that it never comes close to exhaustively exploring the search space in any but the smallest of packer problems, and that the solution returned, while likely to be a *good* solution, is not guaranteed to be the *best* solution.

> **Repeated packer runs are likely to yield a variety of similar solutions near the global optimum; none of these will necessarily *be* the best possible solution.**

## Invoking the packer through the *fixbb* application

To illustrate what the packer does, let's run it.  The simplest way to run the packer is by running the Rosetta *fixbb* application: its sole *raison d'être* is to call the packer on an input pose.  In a terminal, navigate to the demos/public/fixbb directory, and run the following:

```
<path_to_Rosetta_directory>/main/source/bin/fixbb.default.linuxgccrelease -in:file:s 1l2y.pdb -in:file:fullatom -resfile resfile.txt -nstruct 5 >log.txt 2>err.txt &
```

You may need to change "linuxgccrelease", in the above, to whatever is appropriate given your operating system and compiler.

This application packs the side-chains of the input structure (the trp cage mini-protein, 1l2y.pdb).  Five output structures, from five separate runs, are produced.  If you compare these structures to the input structure, you'll find that Rosetta chooses slightly different rotamers for the side-chains, as compared to the input.  This is to be expected, particularly given the discrete nature of rotamers: the truly "best" rotamer might lie between two rotamers tested, and may never be sampled.

Work through the [rest of the demo](../../public/fixbb_design/README).  This teaches about how the packer can be tweaked, both at the commandline and with configuration files called *resfiles*, to control the amount of sampling, the time taken for a run, and the likelihood of converging to the optimal solution.

## The sequence design problem

It is important to note that the packer problem, as described above, makes no assumptions about the nature of the candidate side-chains at each position.  The lists of possibilities can just as easily contain different side-chain identities as it can contain different conformations of the same side-chain.  The packer is therefore a powerful tool (and, indeed, the primary tool in Rosetta) for designing amino acid sequences.

To see the packer design a sequence, open a terminal window, navigate to the demos/public/fixbb_design directory and run the following:

```
<path_to_Rosetta_directory>/main/source/bin/fixbb.default.linuxgccrelease -in:file:s 1l2y.pdb >log.txt &
```

This will produce the output files 1l2y_0001.pdb and score.sc.  If you open 1l2y_0001.pdb in a PDB viewer and compare it to the input file, 1l2y.pdb, you'll see that the sequence has changed considerably.  The rotamers chosen, however, should be interacting with one another reasonably favourably -- that is, there shouldn't be side-chains occupying the same space (clashing), for example.  Note that the only command-line option here is the one to specify our input file; that is, we're not passing any options to the fixbb application in this case to control the behaviour of the packer.  This brings up a very important point:

> **The default behaviour of the packer is to *design* at every position, allowing every rotamer of each of the 20 canonical amino acids.**

If you complete the rest of the [fixbb_design demo](../../public/fixbb_design/README.md), you'll learn how to use resfiles to control the behaviour of the packer, allowing only certain positions to be designed, and only with certain amino acid residue types.

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

There are various setting which can change how the packer behaves, making working with large numbers of rotamers more efficient. (For example, the option `-linmem_ig 10` will change packer settings to be more efficient when using a large number of rotamers per position, though it's less efficient for smaller number of rotamers per position.) However, we can greatly simplify the problem and speed things up by restricting the packer to at most only a few choices at each position with a resfile, as we did in the second part of the second demo:

```
core.pack.pack_rotamers: built 410 rotamers at 16 positions.
```

In general, it is a good idea to limit the packer's options as much as possible in order to keep the search quick and convergent.  This brings up another important concept:

> **Each run of the packer can by controlled with one or more TaskOperations.**

TaskOperations can be passed to the packer in one of several ways, and different TaskOperations modify packer behaviour in different ways.  We've already seen some examples: a user can control packer behaviour with a resfile specified with the "-resfile" option on the commandline, which implicitly invokes the ReadResfile TaskOperation.  Resfiles can control amino acid identity, presence of extra rotamers, and other packer behaviours on a residue-by-residue basis, specified by residue index.  For more information about resfiles and their full features, see the [ReadResfile TaskOperation documentation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ReadResfileOperation) and the documentation for the [resfile syntax](https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/resfiles).  The extra residue options ("-ex1", "-ex2", etc.) also invoke TaskOperations that modify the number of rotamers per amino acid type per position.

One important property of TaskOperations is **commutativity**:

> **TaskOperations can be applied in any order, and still produce the same packer behaviour.**

In order to achieve this, certain TaskOperation-controlled packer behaviours obey *AND* commutativity: if TaskOperation A *and* TaskOperation B *and* TaskOperation C allow the behaviour, then the behaviour will be allowed when A, B, and C are all applied together.  If *any* of A, B, or C prohibits the behaviour, then the behaviour is prohibited when all three are applied.  Allowed canonical amino acid identities at each position obey *AND* commutativity: given TaskOperations A, B, and C, the packer will only design with tyrosine if TaskOperation A *and* B *and* C allow tyrosine at that position.  An easy way to remember this is, "You can only turn canonical residues *off*, and once off, they stay off".  Other TaskOperation functionality obeys *OR* commutativity.  Turning on extra rotamers, for example, occurs if TaskOperation A turns them on *or* TaskOperation B turns them on *or* TaskOperation C turns them on.  (Noncanonical residue identities also obey *OR* commutativity: where design with a particular canonical residue type at a particular position is on by default and can only be turned *off*, design with a particular noncanonical residue type at a particular position is off by default and can only be turned *on*.  Once a noncanonical is on, it stays on.)

As a final note, it's worth being aware that the commutativity of TaskOperations does *not* apply to the MoveMaps that are used to control the minimizer (see the [[minimizer tutorial|minimization]].  In a MoveMap, later commands overrride earlier commands; in a list of TaskOperations, the effects of the TaskOperations combine commutatively.

## Protocols that use the packer

The packer is a fundamental Rosetta algorithm called in the context of many larger protocols.  Many or most protocols call the packer, including *abinitio* (which uses the packer to place and optimize full-atom side-chains after a centroid-mode backbone conformational search) and *relax* (which carries out alternating rounds of packing and energy minimization while ramping the repulsive term in the scoring function).  Nearly every protocol that carries out any sort of design uses the packer for design.

## Calling and controlling the packer from RosettaScripts

In the context of [RosettaScripts](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/RosettaScripts), the packer may be invoked directly using the [PackRotamers](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/PackRotamersMover) mover.  Many other movers, including [FastRelax](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/FastRelaxMover), [FastDesign](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/FastDesignMover), and [Disulfidize](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/DisulfidizeMover) call the packer.  Typically, any Rosetta component that calls the packer can receive one or more TaskOperations to control packer behaviour.  In the RosettaScripts context, [TaskOperations](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/PackRotamersMover) are declared separately in a section of the script preceding the movers that call the packer, and are then passed by name to such movers.

## How the packer algorithm works under the hood (for advanced users)

While a detailed understanding of the workings of the packer is not a strict requirement to use Rosetta, having some idea of what's going on under the hood can help a user to use his or her tools more effectively.  To understand the packer, we must first understand the concept of *Monte Carlo methods*.

Broadly speaking, a Monte Carlo method is any method by which one attempts to solve an optimization problem using random or stochastic moves.  In practice, most Monte Carlo methods in Rosetta make use of the Metropolis-Hastings algorithm, involving iterated steps in which one makes a random move that alters a pose in some way, considers the change in the energy of a pose as a result of the move, and then accepts or rejects that move based on a rule called the Metropolis Criterion.  The Metropolis Criterion states that if the move results in a decrease in energy, it is always accepted, while if it results in an increase in energy, it is accepted with probability equal to e<sup>-ΔE/(k<sub>B</sub>T)</sup>.  Practically, speaking, this means that the larger the increase in energy, the lower the likelihood of accepting the move, while at the same time allowing small increases in energy.

The Metropolis Criterion contains a "temperature factor", k<sub>B</sub>T, that determines how the probability of rejecting a move that increases the energy scales with the magnitude of the increase.  Large values of k<sub>B</sub>T allow moves that result in large increases in energy to be accepted frequently, while small values would only permit moves that result in small increases in energy to be accepted.  While Monte Carlo searches can be carried out with fixed, arbitrarily-selected values of k<sub>B</sub>T, it is often advantageous to vary this value over the course of a simulation, using high values early on to permit the trajectory to climb barriers and escape local energy minima, and then lowering it later in a simulation so that the trajectory "drills down" to the bottom of the lowest-energy well that it has found.  Often, many cycles of ramping k<sub>B</sub>T up and down can allow more complete exploration of the space.  Monte Carlo approaches involving ramping of the k<sub>B</sub>T term are broadly called *simulated annealing methods*.  Such methods often also involve conditional back-stepping to earlier parts of trajectories to escape dead ends, and other adaptive behaviours that increase the probability of converging to the lowest-energy state.

With this understanding of general Monte Carlo and simulated annealing approaches in hand, let's look at what the packer does.  A typical packer run consists of three steps:

1.  TaskOperations are evaluated, and the packer makes a list of possible rotamers at each position.
2.  The packer carries out a precomputation in which all possible pairs of interacting rotamers are enumerated and their pairwise interaction energies are calculated and stored.
3.  The packer carries out a simulated annealing-based search of rotamer combinations, in which moves consist of randomly selecting a position and replacing the current rotamer at that position with a randomly-selected rotamer from the allowed rotamers for that position.  Simulated annealing requires rapid computation of the change in energy resulting from the move.  Because all pairwise interaction energies are precomputed, determining the change in energy of the structure following such a substitution is extremely fast, since it depends only on the internal energies of the old and new rotamers and on their pairwise interaction energies with their neighbours in the pose.  This allows the packer to evaluate hundreds of thousands or millions of moves in seconds, permitting very long trajectories to be carried out very quickly.

There are variants on the above behaviour in which only parts of the interaction network are precomputed and other parts are computed on the fly (what `-linmem_ig` does), but this general scheme is fairly representative.  Note that the above depends heavily on being able to rapidly update the energy as moves are considered:

> **In order to be compatible with the packer, an energy term must either be residue-level pairwise-decomposable, or must otherwise be very fast to compute and update as rotamer substitutions are considered.**
