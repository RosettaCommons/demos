#!/usr/bin/env python2.7
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
#  @author Ray Yu-Ruei Wang, wangyr@u.washington.edu
#
from argparse import ArgumentParser
import os
import shutil
import sys

import denovo_model_building_util

def write_condor_job( args, dir ):

    assert not os.path.exists( args.condor_job_fn )
    job_file = open( args.condor_job_fn, 'w')

    job_file.write("universe      = vanilla\n\n")
    job_file.write("Notify_user   = %s\n" % args.youremail)
    job_file.write("notification  = Error\n\n")
    job_file.write("copy_to_spool = False\n\n")
    job_file.write("Requirements  = ( Memory > 248 )\n")
    job_file.write("Log           = condor.log\n")
    job_file.write("Executable    = /bin/bash\n")
    job_file.write("Output = log\n")
    job_file.write("Error  = err\n\n")
    job_file.write("arguments = ../%s\n\n" % args.run_script )
    job_file.write("Queue 1\n")


def prepare_syd_job( args ):
    submit = open( args.submit_script_fn, "w" )
    submit.write("touch submitted.lock\n")
    num_queues = int( os.popen("grep position %s | wc -l " % args.fragfile ).readline() )

    if args.target_list:
        assert os.path.exists( args.target_list )
        targets = [ int(x.strip()) for x in open( args.target_list, "r" ).readlines() if not x.startswith("#") ]

    curr_dir = os.getcwd()
    for i in range( args.starting_num, args.starting_num + num_queues ):
        # for second run selected residues
        if args.target_list:
            if i not in targets: continue
        sys.stderr.write(".") # just for fun to print out .........

        os.system("mkdir -p %s" % i )
        os.chdir("%s" % i )

        write_condor_job( args, i )

        submit.write("cd ./"+str(i)+"/\n")
        submit.write("condor_submit %s\n" % args.condor_job_fn )
        submit.write("cd ..\n")
        submit.write("sleep 0.1\n")

        os.chdir("%s/" % curr_dir )

    submit.close()
    print "finished at %s/%s" %( i, args.condor_job_fn )
    print "to execute condor_jobs, run 'sh %s'" % args.submit_script_fn



if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument("-t", "--target_list", help="mainly designed for 2nd run when the unassigned regions are defined")
    parser.add_argument("-g", "--fragfile", required=True )
    parser.add_argument("-r", "--run_script", required=True )
    """ for condor_jobs """
    parser.add_argument("--youremail", default="youremail@uw.edu" )
    parser.add_argument("--starting_num", default=1, type=int )
    parser.add_argument("--condor_job_fn", default="condor_job" )
    parser.add_argument("--submit_script_fn", default="submit.sh" )
    args = parser.parse_args()

    assert ("mers" in args.fragfile )
    assert os.path.exists( args.fragfile )

    prepare_syd_job( args )
