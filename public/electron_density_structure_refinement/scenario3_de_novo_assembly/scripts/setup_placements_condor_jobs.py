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
    config = denovo_model_building_util.read_config_file( args.config_file )
    PATH = config.get("path", "demo_dir")
    exe = PATH + "/rosetta/place_fragment_into_density.static.linuxgccrelease"
    db  = PATH + "/rosetta/rosetta_database"

    assert os.path.exists( exe ), exe
    assert os.path.isdir( db ), db

    job_file = open( args.condor_job_fn, 'w')

    job_file.write("universe      = vanilla\n\n")
    job_file.write("Notify_user   = %s\n" % args.youremail)
    job_file.write("notification  = Error\n\n")
    job_file.write("copy_to_spool = False\n\n")
    job_file.write("Requirements  = ( Memory > 248 )\n")
    job_file.write("Log           = condor.log\n")
    job_file.write("Executable    = %s\n\n" % os.path.abspath( exe ))
    job_file.write("Output = log\n")
    job_file.write("Error  = err\n\n")
    job_file.write("arguments = -mute all -pos %s -out:file:scorefile pos_%s_score.sc -out:file:silent pos_%s_silent.out -database %s @../%s/%s \n\n" %( dir, dir, dir, os.path.abspath( db ), args.input_folder, args.flags_fn ))
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


def write_syd_flags( args ):
    buff = open( args.flags_fn, "w" )
    if args.native:
        buff.write("-in:file:native ../%s/%s\n" %( args.input_folder, os.path.basename( args.native )))
    elif args.fasta:
        buff.write("-in:file:fasta ../%s/%s\n" %( args.input_folder, os.path.basename( args.fasta )))
    else:
        stderr.write("ERROR: you need to specify either fasta or native\n")
        exit()

    buff.write("-mapfile  ../%s/%s\n" %( args.input_folder, os.path.basename( args.mapfile )))
    buff.write("-fragfile ../%s/%s\n" %( args.input_folder, os.path.basename( args.fragfile )))
    buff.write("-keep %s\n" % args.poses_to_keep )
    buff.write("-cycles 3\n")
    buff.write("-min_pack_min\n")
    buff.write("-skip_edge 1\n")
    buff.write("-no_density_score %s\n" % args.no_density_score )
    buff.write("-chemical:exclude_patches LowerDNA UpperDNA Cterm_amidation SpecialRotamer protein_cutpoint_upper protein_cutpoint_lower VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm\n")
    buff.close()


def prepare_input_folder( args ):
    write_syd_flags( args )
    assert not os.path.exists( args.input_folder )

    os.mkdir( args.input_folder )

    if args.native:
        shutil.copy( args.native, args.input_folder )
    else:
        shutil.copy( args.fasta, args.input_folder )

    shutil.copy( args.mapfile, args.input_folder )
    shutil.copy( args.fragfile, args.input_folder )
    shutil.move( args.flags_fn, args.input_folder )


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('--config_file', default="../denovo_model_building_scripts.cfg", help='')
    parser.add_argument("-t", "--target_list", help="mainly designed for 2nd run when the unassigned regions are defined")
    parser.add_argument("-n", "--native", help="")
    parser.add_argument("-f", "--fasta", help="")
    parser.add_argument("-m", "--mapfile", required=True)
    parser.add_argument("-g", "--fragfile", required=True )
    parser.add_argument("--poses_to_keep", default=2000 )
    parser.add_argument("--no_density_score", default=3.0 )
    parser.add_argument("--flags_fn", default="placements_flags" )
    parser.add_argument("--input_folder", default="input_files" )
    """ for condor_jobs """
    parser.add_argument("--youremail", default="youremail@uw.edu" )
    parser.add_argument("--starting_num", default=1, type=int )
    parser.add_argument("--condor_job_fn", default="condor_job" )
    parser.add_argument("--submit_script_fn", default="submit.sh" )
    args = parser.parse_args()

    assert ("mrc" in args.mapfile) and ("mers" in args.fragfile )
    assert os.path.exists( args.mapfile) and os.path.exists( args.fragfile )

    if args.native:
        assert ("pdb" in args.native ) and os.path.exists( args.native )
    elif args.fasta:
        assert ("fasta" in args.fasta ) and os.path.exists( args.fasta )
    else:
        stderr.write("ERROR: you need to specify either fasta or native\n")
        exit()


    prepare_input_folder( args )
    prepare_syd_job( args )
