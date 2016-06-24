# Multistate Design of Antibodies

KEYWORDS: DESIGN GENERAL

This demo includes input files for running multistate design using 10 states.
It can either be run on a single processor using the
```
   mpi_msd.default.linuxgccrelease (or, more generally the mpi_msd.default.{os}{compiler}{release/debug} executable)
```
   using the provided input, you can test this:
   
   (where `$ROSETTA3`=path-to-Rosetta/main/source

```
   $> $ROSETTA3/bin/mpi_msd.linuxgccrelease @flags
```

or by distributing those states over multiple processors using the
```
   mpi_msd.mpi.linuxgccrelease (or, more generally the mpi_msd.mpi.{os}{compiler}{release/debug} executable)
```
and by launching the job with between 1 and 10 CPUs.  Using 5 CPUs would place two states
on each processor.

Files to look at first:
```
1USM_het.flags:   contains the set of flags read by the command line in the command/command_mpi files.
fitness.daf:      the file which declares the set of states included in this design task, and the fitness function itself
entity.resfile:   the resfile which declares the accessible regions of sequence space
```
