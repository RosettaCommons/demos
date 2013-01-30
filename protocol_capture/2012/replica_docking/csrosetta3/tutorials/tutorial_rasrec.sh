#!/bin/bash
###
###
### This file is part of the CS-Rosetta Toolbox and is made available under
### GNU General Public License
### Copyright (C) 2011-2012 Oliver Lange
### email: oliver.lange@tum.de
### web: www.csrosetta.org
###
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program.  If not, see <http://www.gnu.org/licenses/>.
###
###


# ----------------------------------------------------------------------------------------------------------------
# PREAMBLE
# the code in the PREAMBLE is specific to make these scripts work as tutorials.
# this code is not usually carried out by a user who performs the steps explained in the TUTORIAL part of this script
#
#
function tail_log() {
#parameter file, N
echo tail of output: $1
echo + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
tail -n $2 $LOGS/$1
echo + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
echo
}

# overwrite to keep tutorial targets out of users target-lib,
# note that a relative path here would mean that setup_target and setup_run
# have to be started from the same directory. With absolute paths this restriction does not apply
CS3_BENCH_TARGETLIB=$(pwd)/tutorial_targets


if [ ! -e inputs ]; then
    echo
    echo ERROR: the inputs directory of the CSROSETTA toolbox is required
    echo Solution1: either start the tutorials in the folder they are distributed in
    echo Solution2: issue the command 'ln -s <patch_to_csrosetta3>/tutorials/inputs' in your local folder
    echo Trying second solution automatically: ln -s $( dirname $( dirname `which setup_target` ) )/tutorials/inputs
    echo
    ln -fs $( dirname $( dirname `which setup_target` ))/tutorials/inputs
    if [ -e inputs ]; then
	if [ -e inputs/2jrm_trim.fasta ]; then
	    echo 'Succeeded in creating a link. Starting the tutorial'
	    echo
	else
	    echo 'Failed in creating the link'
	    echo
	    exit 1
	fi
    else
	echo 'Failed in creating the link'
	echo
	exit 1
    fi
fi

#overwrite to keep tutorial targets out of users target-lib
CS3_BENCH_TARGETLIB=$(pwd)/tutorial_targets

mkdir -p run_rasrec_logs
LOGS=$(pwd)/run_rasrec_logs

# END_PREAMBLE
# -----------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------
# START TUTORIAL
#
# NOTE 1: The redirecting of output with > log file && tail_log || exit
#    is also tutorial specific as a user you probably run these scripts interactively
#
echo
echo TUTOR: Tutorial \'rasrec\' is starting...
echo

#1 prepare target from scratch
echo TUTOR: Setup Target with fragments, fasta and chemical shifts...
echo
echo CMD: setup_target -target 2jrm_trim -method rasrec -frags inputs/2jrm_trim.frags*dat.gz \
    -fasta inputs/2jrm_trim.fasta -cs inputs/2jrm_trim.tab
echo Redirecting output to run_rasrec_logs/setup_target.1.log
setup_target -target 2jrm_trim_2 -method rasrec -frags inputs/2jrm_trim.frags*dat.gz -fasta \
             inputs/2jrm_trim.fasta -cs inputs/2jrm_trim.tab > $LOGS/setup_target.1.log \
&& tail_log setup_target.1.log 15 || exit 1

#2 change method on existing target
echo TUTOR: Setup Target by transferring from an existing Setup
echo
echo CMD: setup_target -target 2jrm_trim -transfer_method abrelax -method rasrec
echo Redirecting output to run_rasrec_logs/setup_target.2.log
setup_target -target 2jrm_trim -transfer_method abrelax -method rasrec > $LOGS/setup_target.2.log \
&& tail_log setup_target.2.log 15 || exit 1

#3. setup a ROSETTA run
if [ -e run_rasrec ]; then
    echo
    echo RUN directory run_rasrec exists already... Just checking for completed runs
    echo
    CONTINUE=yes
else
    CONTINUE=no
fi

if [ "$CONTINUE" == "no" ]; then
 echo TUTOR: Generate a CS-ROSETTA RUN directory
 echo
 echo CMD: setup_run -target 2jrm_trim -method rasrec -dir run_rasrec [-job slurm/moab]
 echo Redirecting output to run_rasrec_logs/setup_run.log
 setup_run -target 2jrm_trim -method rasrec -dir run_rasrec > $LOGS/setup_run.log && \
 tail_log setup_run.log 5 || exit 1

#NOTE: the cycle-factor is set here very low to have the tutorial finish fast.
#NOTE: Default setting for RASREC is 2
fi

#select here to run test or production (test is much faster, but doesn't produce meaningful results)
RDIR=test
RUN=test
# RUN=production
# RDIR=run
echo INFO: Selected to run $RUN version of protocol -- uncomment/comment in script to select between test and production

#change to run directory
cd run_rasrec/2jrm_trim/$RDIR

if [ "$CONTINUE" == "no" ]; then
 # start a 48-core job (200-500 cores are the recommended range)
 echo
 echo TUTOR: 'Start a multiprocessor (48 core) job using the SLURM queuing system ...'
 echo CMD: sbatch -n 48 $RUN.slurm.job
 type sbatch >/dev/null 2>&1  && \
 sbatch -n 48 $RUN.slurm.job && queue=slurm || \
     echo 'sbatch not found, you probably have a different queuing system' && queue=none
 echo
 if [ "$queue" == "none" ]; then
   echo TUTOR: 'Start a multiprocessor (64 core) job using the MOAB (as e.g., on JUROPA/Juelich) queuing system ...'
   echo CMD: msub -l nodes=4:ppn=16 $RUN.slurm.job
   type msub >/dev/null 2>&1  && \
   msub -l nodes=4:ppn=16 $RUN.slurm.job && queue=moab || \
   echo 'msub not found, you probably have a different queuing system'
 fi
fi

if [ "$queue" == "slurm" ]; then
    echo TUTOR: check output in slurm-output file
    tail -n 5 slurm*out
fi

echo
echo TUTOR: check log-files in directory logs_'<jobid>/'
echo TUTOR: the most important log is logs_'<jobid>/log_2', which contains the Pool-Manager process
echo TUTOR: log_0 is from the dedicated IO process and usually contains no errors
echo TUTOR: log_1 is from the job-distributor process and also usually without errors
echo TUTOR: log_3 and higher is from the worker processes. These are all identical.
echo TUTOR: look for errors in log_2 and log_3. If you get a specific message that job crashed with rank X,
echo TUTOR: look in log_X.
echo

echo TUTOR: check for completion - there should be a directory fullatom_pool_stage8
if [ -e fullatom_pool_stage8 ] && [ -e fullatom_pool_stage8/decoys.out ]; then
    echo $RUN 'run is completed'
else
    echo $RUN 'run is not completed'
    echo 'Wait (several hours/days)'
    if [ -e batch_000001 ]; then
	echo Batches generated so far:
	ls -d batch*
    fi
    if [ -e centroid_pool_stage1 ]; then
	echo Stages finished
	ls -d *pool_stag*
    fi
    exit 0
fi

#final decoys are in subdirectory fullatom_pool
cd fullatom_pool

#extract rms and score for ploting
echo TUTOR: Extract scores and low-energy decoys
echo
echo CMD: 'extract_scores decoys.out score rms chem_shifts description > rms_score.txt'
extract_scores decoys.out rms score chem_shifts description > rms_score.txt
echo this created the file rms_score.txt that has this content:
cat rms_score.txt

#extract low-scoring decoys '(usually 10 of 20.000-50.000)'
echo
echo CMD: 'extract_decoys decoys.out -score 10 > low_10.out'
extract_decoys decoys.out -score 10 > low_10.out
echo this created the file $( ls low_10.out )

#convert to PDB format
echo
echo TUTOR: Convert from silent file to PDB
echo CMD: make multi-model PDB file from 10 lowest energy structures
pack_pdbs -silent low_10.out > final.pdb

echo
echo TUTOR: get individual files name as in low_10.out
echo 'CMD: cat final.pdb | unpack_pdbs -remark ROSETTA-TAG'
cat final.pdb | unpack_pdbs -remark ROSETTA-TAG
echo this created files: $( ls bat*pdb )
echo
echo TUTOR: get individual files named model_01...model_10
echo 'CMD: cat final.pdb | unpack_pdbs'
cat final.pdb | unpack_pdbs
echo this created files: $( ls model*pdb )
echo
echo These files can be found in  $( pwd )
cd ../../../..

