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
#how fast can you read?
SLEEP='sleep 5'

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

function cmd() {
#parameters cmd, log, N, exit_cmd
log=$2
N=$3
if [ $# -lt 4 ]; then
  exit_cmd='exit 1'
else
  exit_cmd="$4"
fi

echo CMD: $1
echo Redirecting output to $( echo $LOGS | sed s@$(pwd)/@@ )/$log.log
$1 > $LOGS/$log.log \
&& tail_log $log.log $N || echo $exit_cmd
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

mkdir -p run_restraints_logs
LOGS=$(pwd)/run_restraints_logs

# END_PREAMBLE
# -----------------------------------------------------------------------------------------------------


# ------------------------------------------------

# -----------------------------------------------------------------------------------------------------
# Start TUTORIAL
#
# NOTE 1: The redirecting of output with > log file && tail_log || exit
#    is also tutorial specific as a user you probably run these scripts interactively
#
echo
echo TUTOR: Tutorial \'restraints\' is starting...
echo

#1 prepare target from scratch
echo TUTOR: Setup Target with fragments, fasta and chemical shifts...
echo
cmd 'setup_target -target SgR145 -method rasrec -frags inputs/SgR145.frags*dat.gz -fasta inputs/SgR145.fasta -cs inputs/SgR145.tab' setup_target.basic 15


#2 add RDC data
echo TUTOR: Setup Target by transferring from an existing Setup
cmd 'setup_target -target SgR145 -method rasrec -rdc inputs/SgR145*.rdc' setup_target.rdc 15
$SLEEP

#3 add single medium RDC data
echo TUTOR: Setup a labelled target with the first RDC medium
cmd 'setup_target -target SgR145 -method rasrec -transfer_label standard -label rdc_med1 -rdc inputs/SgR145.med1.rdc' setup_target.rdc1 15

echo TUTOR: Setup a labelled target with the other RDC medium
cmd 'setup_target -target SgR145 -method rasrec -transfer_label standard -label rdc_med2 -rdc inputs/SgR145.med2.rdc' setup_target.rdc2 15
echo TUTOR: Note that the data files SgR145.med1.rdc and SgR145.med2.rdc are uploaded only once. The md5 checksum shows that the content is
echo unaltered the second time the file is uploaded, otherwise setup_target would complain
$SLEEP

cp inputs/SgR145.med2.rdc .
echo 'CHANGED LOCAL COPY' >> SgR145.med2.rdc
echo TUTOR: Attempting to upload an altered version of the med2 file with the same name:
cmd 'setup_target -target SgR145 -method rasrec -transfer_label standard -label rdc_med2 -rdc SgR145.med2.rdc' setup_target.rdc2 15 'echo "TUTOR: As you see above, this caused an error"'
echo waiting for you to read...
$SLEEP

echo thanks.

echo TUTOR: Setup a run from '(NOE) distance restraints given in CYANA upl-format:'
cmd 'setup_target -target SgR145 -method rasrec -cyana_upl inputs/SgR145.final.upl inputs/SgR145.manual.upl' setup_target.upl 15
$SLEEP

rm -Rf run_restraints

echo TUTOR: When generating the run these restraints in the CYANA upl file
echo TUTOR: 'are separated into high-quality (noQF) and low-quality (QFall),'
echo TUTOR: timmed to match fasta sequence and converted to Rosetta format
cmd 'setup_run -target SgR145 -method rasrec -dir run_restraints' setup_run.upl 3
echo TUTOR: The trimming can be seen in the output of the setup_run command:
echo '======================== somewhere in setup_run.upl.log ============================='
head -n 70 $LOGS/setup_run.upl.log | grep -v subst | grep -v inputs | tail -n 20
echo '====================================================================================='
echo
echo TUTOR: the generated files can be found in run_restraints/SgR145/inputs/nmr_data
echo LIST DIR: $( ls run_restraints/SgR145/inputs/nmr_data)
echo
echo TUTOR: and are included in the simulation via ConstraintsClaimer in the TopologyBroker File:
echo SHOW cat run_restraints/SgR145/inputs/setup_init.tpb:
echo
cat run_restraints/SgR145/inputs/setup_init.tpb
$SLEEP

echo TUTOR: 'For each restraint input file a ConstraintClaimer BLOCK is present.'
echo TUTOR: 'the entries in these blocks from CLAIMER...END_CLAIMER regulate the exact use of the restraints'
echo TUTOR: 'FULLATOM / NO_CENTROID : restraints are only used for full-atom models'
echo TUTOR: 'CENTROID : restraints are only used for centroid models (no side-chain atoms but the CB)'
echo TUTOR: 'COMBINE_RATIO 2 : randomly paired restraints are combined into ambiguous restraints, this drastically'
echo TUTOR: '   reduces impact of wrong restraints'
echo TUTOR: 'FILTER_WEIGHT: how much does the individual set of restraints contribute to the total energy when structures'
echo TUTOR: '   are selected into the pool (global weight is in run/nmr_pool_patch)'
echo TUTOR: 'FILTER_NAME: name used for score-column in SCORE-line of output decoys'



