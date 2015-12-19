#!/usr/bin/python
from os import system
from os.path import exists
from os.path import abspath
import string

if( exists("master_log.out") ): system("rm master_log.out")
if( exists("master_log.err") ): system("rm master_log.err")
if( exists("SLAVE_JOBS/") ): system("rm -rf SLAVE_JOBS/")

master_tag = abspath( "MASTER" ).replace("/","_")

command_act="bsub -W 144:0 -o master_log.out -e master_log.err " + (" -J %s " %(master_tag)) + " ~/SWA_RNA_python/SWA_dagman_python/dagman/DAG_continuous.py  -j 500 rna_build.dag"
system(command_act)
