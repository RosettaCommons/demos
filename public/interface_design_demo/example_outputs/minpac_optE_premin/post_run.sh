#!/usr/bin/env bash
set -e

# stuff you might want to change:

# number of trajectories to use for sequence recovery
n_traj=5

# list of native pdbs for comparison
# nat_pdb_list path is relative to the output directory
# must have exactly one pdb in $nat_pdb_list for each one given for inputs
nat_pdb_list=../../../../inputs/min_nats/nats.list

benchamre_file=benchmark.txt

# end of stuff you might want to change


# The $jobname assignment line is set by run.sh. Putting an exit command in the
# default line ensures run.sh has executed and overwritten the line.
# original: jobname=; echo "run.sh needs to be executed before post_run.sh"; exit
jobname=minpac_optE_premin

# infpro.sh needs to be run with the output directory as working directory
# (It's a modified version of a script in my $PATH)
cd $jobname

../../../analyze/infpro.sh $n_traj $nat_pdb_list
# The script could also have been run without arguments, it has default
# values equal to those passed into it above
# ../../../analyze/infpro.sh

# end of standard run analysis


# commands below generate a file, benchmark.txt, that summarizes performance
# metrics of the run

# Generate a benchmark.txt file
CPU_secs=$(grep -m 1 'CPU time' ${jobname}* | sed 's/^.*CPU time[^0-9]*\([0-9]\+\(\.[0-9]\+\)\) sec.*$/\1/g')
CPU_hours=$(echo "scale=3; ${CPU_secs}/3600" | bc)
tmp=$(grep -m 1 'Total' sequencerecovery.txt)
seq_rec_core=$(echo $tmp | awk '{ print $4}')  # I'm unaware of a good way to
seq_rec_tot=$(echo $tmp | cut -d ' ' -f 7)     # extract fields from strings,
seq_rec_surf=$(echo $tmp | awk '{ print $11}') # here are 2 ugly ones

# Grab some of the run parameters
num_prots=$(wc -l start_struct.list | awk '{print $1}')
weights=$(grep -m 1 '\-score:weights' ${jobname}.run* | sed 's/^.*-score:weights \([[:graph:]]\+\) .*$/\1/g')
nstructs=$(grep -m 1 '\-out:nstruct [0-9]' ${jobname}.run* | sed 's/^.*-out:nstruct \([0-9]\+\).*$/\1/g')

cd ..
echo "Benchmark file" > benchmark.txt
echo "at $(date)" >> benchmark.txt
echo -e "in $(pwd)\n" >> benchmark.txt
echo "Weights used: ${weights}" >> benchmark.txt
echo "Number of input structures: ${num_prots}" >> benchmark.txt
echo "Number of design trajectories: ${nstructs}" >> benchmark.txt
echo "CPU hours: ${CPU_hours}" >> benchmark.txt
echo "Overall interface sequence recovery: ${seq_rec_tot}" >> benchmark.txt
echo "Core interface sequence recovery: ${seq_rec_core}" >> benchmark.txt
echo "Surface interface sequence recovery: ${seq_rec_surf}" >> benchmark.txt
