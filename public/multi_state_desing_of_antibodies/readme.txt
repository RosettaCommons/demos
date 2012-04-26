This demo includes input files for running multistate design using 10 states.
It can either be run on a single processor using the

   mpi_msd.default.linuxgccrelease (or, more generally the mpi_msd.default.{os}{compiler}{release/debug} executable)

or by distributing those states over multiple processors using the

   mpi_msd.mpi.linuxgccrelease (or, more generally the mpi_msd.mpi.{os}{compiler}{release/debug} executable)

and by launching the job with between 1 and 10 CPUs.  Using 5 CPUs would place two states
on each processor.

Files to look at first:

command :         contains an example command line which you will have to modify for your system
command_mpi:      contains an example command line for running the mpi executable which you will have to modify for your system
1USM_het.flags:   contains the set of flags read by the command line in the command/command_mpi files.
fitness.daf:      the file which declares the set of states included in this design task, and the fitness function itself
entity.resfile:   the resfile which declares the accessible regions of sequence space
