Heterodimeric Antibody Design using Multistate Design
=====================================================

KEYWORDS: ANTIBODIES DESIGN

Author: Andrew Leaver-Fay

---

Multistate design considers the impact that a sequence has on multiple structures
(states) simultaneously to rule one sequence more favorable (fit) for a particular
purpose than another sequence.  In the case of this protocol, I'm attempting to
design a heterodimer starting with a homodimer.  Multistate design needs to favor
the binding interactions for the heterodimeric species (AB) while disfavoring
the binding interactions for the homodimeric species (AA and BB).  This objective
is encoded in an input-file fitness function.

A new multistate protein design protocol has been implemented in Rosetta.  In an outer loop,
a genetic algorithm explores sequence space; in an inner loop, fixed-sequence rotamer repacking
finds a low-energy rotamer assignments for each of arbitrarily many fixed-backbone states.
(Fixed sequence repacking can be distributed across multiple CPUs with MPI).  The energies for
the set of states are used to define a fitness for a sequence; the function that converts from
the state energies to a fitness is definable in an input file and may be arbitrarily complex.
This protocol has been applied toward the design of bispecific antibodies by redesigning the Fc
interface; the protocol sought sequences favoring the AB heterodimer while disfavoring the AA
and BB homodimers. 

Running the protocol
--------------------

#### Essential Flags

    -entity_resfile  <fname>
    -fitness_file    <fname>
    -ms::generations <int>
    -ms::pop_size    <int>
    -ms::fraction_by_recombination <float>

The entity resfile defines the number of entity elements in the design task, and defines the
amino acid and rotamer search space for each entity element.  The entity resfile file format
is simply a resfile that's proceeded by one line containing one integer, the number of elements
in the entity.

The fitness file defines the set of states which are to be optimized, and the fitness function
that determines the fitness for an entity, given the energies of each of the states after
they've been repacked using the sequence encoded in that entity.  There are six available
commands in the fitness file; STATE, STATE_VECTOR, VECTOR_VARIABLE, SCALAR_EXPRESSION,
VECTOR_EXPRESSION, ENTITY_FUNCTION, and FITNESS.  Each command must be on its own line.
The syntax for command X may be found on in the @details section preceeding the
function definition for:

    void DynamicAggregateFunction::process_X_command

in the file

    /mini/src/devel/pack_daemon/DynamicAggregateFunction.cc

#### Important Flags

    -ms::generations <int>

I have found that the number of generations should be ~ 15 x N where N
is the number of entity elements. In this example, I have 7 residues on each side of the
interface being designed, so I have 14 entity elements, and therefore run for 210 generations.

    -ms::pop_size <int>

I have found that, using ms::generations <15\*N>, the population size of 100 is good.

    -ms::fraction_by_recombination <float>

The fraction of crossover events; the fraction mutated by 
point mutations is 1 - this number.  I have found that high point mutation rates are preferable, and
typically set the fraction by recombination to 5% or lower.

#### Less commonly used flags

    -ms::numresults <int>

The number of results to output, sorted by increasing (worsening) fitness.
Typically there are many sequences near the top sequence that are not very different in sequence space
nor in fitness.  The default is to output only the pdbs relating to the entity with the best fitness.

All of the flags that control the initialization of a packer-task may be used, but their use is
discouraged in favor of specifying behavior for residues in either the entity-resfile or the secondary
resfiles.

#### Example Rosetta Command Line:

    path/to/mini/bin/mpi_msd.linuxgccrelease -entity_resfile input_files/entity.resfile -fitness_file input_files/fitness.daf  -ms::pop_size 100 -ms::generations 210  -ms::numresults 1 -no_his_his_pairE -ms::fraction_by_recombination 0.02 -database /path/to/minirosetta_database

#### Example Overall Command Line (if overall protocol is run via a script or other program):

One MSD run is insufficient for multiple reasons: 1) one design run is never enough, 2) your fitness function
probably has one or more parameters that need sweeping, and should be swept through, 3) if you're modeling negative
states, then you probably need to iteratively generate those negative states.  You probably also
need scripts to control the submission of the MSD jobs to an MPI cluster.  That will vary from cluster to cluster.

In the case of the heterodimer design, I also used two more rosetta applications:

    docking_protocol.linuxgccrelease

to re-dock the homodimeric species following their output from the multistate design protocol
with the flags

    -s <pdbname>
    -database /path/to/minirosetta_database
    -nstruct 20
    -docking:docking_local_refine 1
    -no_his_his_pairE


and

    InterfaceAnalyzer.linuxgccrelease

with the flags

    -s <pdbname>
    -database /path/to/minirosetta_database
    -jd2::no_output
    -jumpnum 1
    -overwrite
    -is_compute_hbond_unsat true
    -is_compute_packstat true
    -mute protocols.toolbox
    -no_his_his_pairE

to measure the interface energy, the buried surface area (dSASA), the "binding energy density" (interface energy / dSASA),
and the number of buried unsatisfied hydrogen bonds.

Versions
--------

* Committed to trunk at revision 36337
* InterfaceAnalyzer was used with version "32988:33373M"
* docking_protocol was used with version "34393:34892M"

Other Comments
--------------
Every state which contributed in some way to the fitness for a particular entity is output
at the conclusion of the MSD run. The state is output with the rotamer assignment computed
when the entity's fitness was first evaluated (states are not repacked prior to being output,
so if you think there is something fishy in a rotamer packing, you can look at the rotamer
assignment). For example, if you used the "vmin" function to select the state 
from a state vector with the lowest energy, then only the state with the lowest 
energy is output.

Output pdbs are named with "msd\_output\_", the rank of the source entity 
(1..numresults), the name for the state variable or the state-vector variable 
(this name comes from the fitness file) and a ".pdb".

