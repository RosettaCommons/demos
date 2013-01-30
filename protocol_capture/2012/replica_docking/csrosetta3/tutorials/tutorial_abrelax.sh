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

mkdir -p run_abrelax_logs
LOGS=$(pwd)/run_abrelax_logs

# END_PREAMBLE
# -----------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------
# START TUTORIAL
#
# NOTE 1: the way I switch between directories inputs and . is not necessarily the best practice but
# A) serves to demonstrate that you can setup the same target from different locations
# B) allows input to remain a readonly directory so that you can run tutorials locally even when the CSROSETTA toolbox has been installed
#    on your system by another user.
# NOTE 2: The redirecting of output with > log file && tail_log || exit
#    is also tutorial specific as a user you probably run these scripts interactively
#
echo
echo TUTOR: Tutorial \'abrelax\' is starting...
echo

#1. prepare target
cd inputs
echo TUTOR: Setup Target with fragments, fasta and chemical shifts...
echo
echo CMD: setup_target -target 2jrm_trim -method abrelax -frags 2jrm_trim.frags*dat.gz -fasta 2jrm_trim.fasta -cs 2jrm_trim.tab
echo Redirecting output to run_abrelax_logs/setup_target.1.log
setup_target -target 2jrm_trim -method abrelax -frags 2jrm_trim.frags*dat.gz -fasta 2jrm_trim.fasta -cs 2jrm_trim.tab > $LOGS/setup_target.1.log && \
tail_log setup_target.1.log 15 || exit 1

cd ..

#2. optional step, add native PDB for RMSD computation
echo TUTOR: Add a native PDB for on-the-fly RMSD compuation
echo
echo CMD 'renumber_pdb inputs/2jrmA.pdb -fasta inputs/2jrm_trim.fasta > 2jrm_trim.pdb'
renumber_pdb inputs/2jrmA.pdb -fasta inputs/2jrm_trim.fasta > 2jrm_trim.pdb
echo
echo CMD setup_target -target 2jrm_trim -method abrelax -native 2jrm_trim.pdb
echo Redirecting output to run_abrelax_logs/setup_target.2.log
setup_target -target 2jrm_trim -method abrelax -native 2jrm_trim.pdb > $LOGS/setup_target.2.log && \
tail_log setup_target.2.log 15 || exit 1

rm 2jrm_trim.pdb

#3. setup a ROSETTA run
if [ -e run_abrelax ]; then
    echo
    echo RUN directory run_abrelax exists already... Delete and start over
    echo DO: rm -Rf run_abrelax
    echo
    exit 1
fi
echo TUTOR: Generate a CS-ROSETTA RUN directory
echo
echo CMD: setup_run -target 2jrm_trim -method abrelax -dir run_abrelax -cycle_factor .01 -nstruct 10 -job interactive -relax
echo Redirecting output to run_abrelax_logs/setup_run.log
setup_run -target 2jrm_trim -method abrelax -dir run_abrelax -cycle_factor .01 -nstruct 20 -job interactive -relax > $LOGS/setup_run.log && \
tail_log setup_run.log 5 || exit 1

#NOTE: the cycle-factor is set here very low to have the tutorial finish fast.
#NOTE: Default setting for RASREC is 2 and
#NOTE: for abrelax the default setting is 10

#change to run directory
cd run_abrelax/2jrm_trim/run

## single processor run
# source production.interactive.job

# or using multiple cores
echo TUTOR: 'Start a multiprocessor (22 core) job ...'

echo CMD: source production.interactive.job -n 22
echo Redirecting output to abrelax.log
source production.interactive.job -n 22 > $LOGS/abrelax.log && \
tail_log abrelax.log 15 || exit 1


#extract rms and score for ploting
echo TUTOR: Extract scores and low-energy decoys
echo CMD: 'extract_scores decoys.out score rms chem_shifts description > rms_score.txt'
extract_scores decoys.out rms score chem_shifts description > rms_score.txt
echo this created the file rms_score.txt that has this content:
cat rms_score.txt


#extract low-scoring decoys '(usually 10 of 20.000-50.000)'
echo CMD: 'extract_decoys decoys.out -score 10 > low_10.out'
extract_decoys decoys.out -score 10 > low_10.out

#convert to PDB format
echo
echo TUTOR: Convert from silent file to PDB
echo CMD: make multi-model PDB file from 10 lowest energy structures
pack_pdbs -silent low_10.out > final.pdb

echo
echo TUTOR: get individual files name as in low_10.out
echo 'CMD: cat final.pdb | unpack_pdbs -remark ROSETTA-TAG'
cat final.pdb | unpack_pdbs -remark ROSETTA-TAG
echo this created files: $( ls S*pdb )
echo
echo TUTOR: get individual files named model_01...model_10
echo 'CMD: cat final.pdb | unpack_pdbs'
cat final.pdb | unpack_pdbs
echo this created files: $( ls model*pdb )
echo
echo These files can be found in  $( pwd )
cd ../../..


