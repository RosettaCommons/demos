#!/bin/bash

#PBS -N biox3_scratch_users_rhiju_projects_stepwise_test_EETI-II_2it7_new_test_SLAVE_JOBS_99

#PBS -o SLAVE_JOBS/99/slave_jobs.out_QSUB

#PBS -e SLAVE_JOBS/99/slave_jobs.err_QSUB

#PBS -l walltime=48:00:00

cd $PBS_O_WORKDIR

cat $PBS_NODEFILE > SLAVE_JOBS/99/nodefile.txt

/home/rhiju/SWA_dagman_python2/SWA_dagman_slave.py SLAVE_JOBS/99 >SLAVE_JOBS/99/slave_jobs.out 2>SLAVE_JOBS/99/slave_jobs.err

