#!/usr/bin/tcsh

c4 -cvs warhol-01..15 ls /vol/ek/ravehb
echo RUN FLAGS:
cat output_files/flags
mpirun  -np 112 -machinefile warhol.mpi --mca plm_rsh_agent "rsh : ssh -X -n"  \
  tcsh ./run_example
