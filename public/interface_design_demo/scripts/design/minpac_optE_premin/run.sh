#!/usr/bin/env bash
set -e

jobname="minpac_optE_premin"

# Set the job name in the post_run.sh script.
sed -i "s/^jobname=.*$/jobname=${jobname}/g" post_run.sh

# Run the script that actually executes the binary. Pass it the job name, the
# list of input pdbs, and the Rosetta script.

runner_with_options=${jobname}.sh
input_pdb_list=selected_chains.list
rosetta_script=${jobname}.xml

./$runner_with_options $jobname $input_pdb_list $rosetta_script
