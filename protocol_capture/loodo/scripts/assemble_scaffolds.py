#!/usr/bin/python


####################################################################################
# For assembling CapHits and the two halves of the bottom domain into the same PDB file
####################################################################################


from rosetta import *
import os, sys, getopt

def main(argv):
    
    loodo_bot = ''
    loodo_ins_begin = ''
    caphits_listfile = ''
    caphits_directory = ''
    init_flags = ''
    
    try:
        opts, args = getopt.getopt(argv, "hb:i:l:d:f:", ["bot_file=", "ins_begin=", "list_of_caphits=", "caphits_directory", "init_flags"])
    except getopt.GetoptError:
        print('assemble_scaffolds.py -b <bottom_domain_pdbfile> -i <insertion_site_beginning_residue> -l <caphit_RTs_listfile> -d <directory_containing_caphits> -f <PyRosetta_init_flags>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print 'assemble_scaffolds.py -b <bottom_domain_pdbfile> -i <insertion_site_beginning_residue> <caphit_RTs_listfile> -d <caphits_directory_OPTIONAL> -f <PyRosetta_init_flags>'
            sys.exit()
        elif opt in ("-b", "--bot_file"):
            loodo_bot = arg
        elif opt in ("-i", "--ins_begin"):
            loodo_ins_begin = arg
        elif opt in ("-l", "--list_of_caphits"):
            caphits_listfile = arg
        elif opt in ("-d", "--caphits_directory"):
            caphits_directory = arg
        elif opt in ("-f", "--init_flags"):
            init_flags += '-'
            init_flags += arg
            init_flags += ' '

    if loodo_bot=='' or loodo_ins_begin=='' or caphits_listfile=='':
        sys.exit("Missing one or more required inputs. Use assemble_scaffolds.py -h to see usage information.")

    print("Arguments:\n \
          Parent Domain (Bot): {bot}\n \
          Insertion Site (Beginning Residue): {ins_begin}\n \
          CapHits List File: {caphitsfile}\n \
          CapHits Directory: {caphitsdir}\n \
          PyRosetta Init Flags: {init_flags}".format( bot=loodo_bot, 
                                                     ins_begin=loodo_ins_begin, 
                                                     caphitsfile=caphits_listfile, 
                                                     caphitsdir=caphits_directory,
                                                     init_flags=init_flags))
    init(init_flags)

    with open('{}'.format(caphits_listfile)) as caphitfile:  ##CapHit_RT.txt
        caphitlines = caphitfile.readlines()

    bot = pose_from_pdb(loodo_bot)

    ### Assemble them with bottom domain pieces.
    print("Number of Scaffolds: {}".format(len(caphitlines)))
    
    for cap in caphitlines:
        
        capfile = caphits_directory + '/' + cap.split()[0]
        base_name = capfile.split('/')[-1].replace('.pdb','')
        lAB = base_name.split('_')[1]

        link_dir = ''
        if caphits_directory != '':
            link_dir = caphits_directory
        
        linkA = "{link_dir}/LinkerA_{lAB}_{linkerAinfo}.pdb".format( link_dir=link_dir, lAB=lAB, linkerAinfo='_'.join(base_name.split('_')[5:]))
        linkB = "{link_dir}/LinkerB_{lAB}_{linkerBinfo}.pdb".format( link_dir=link_dir, lAB=lAB, linkerBinfo='_'.join(base_name.split('_')[2:4]))
        
        print(capfile)

        cap_pose = pose_from_pdb(capfile)
        linkerA_pose = pose_from_pdb(linkA)
        linkerB_pose = pose_from_pdb(linkB)
        
        #zincs = []

        assembled = Pose()
        assembled.append_residue_by_jump(bot.residue(1),0,"","",0)
        #print "Bot1"
        for i in range(2, int(loodo_ins_begin)-1):
            assembled.append_polymer_residue_after_seqpos(bot.residue(i), assembled.total_residue(),0)
        #print assembled.sequence()
        #print "LinkerA"
        for i in range(1, linkerA_pose.total_residue()+1):
            assembled.append_polymer_residue_after_seqpos(linkerA_pose.residue(i), assembled.total_residue(),0)
        #print assembled.sequence()
        #print "Cap"
        for i in range(3, cap_pose.total_residue()-1):
            assembled.append_polymer_residue_after_seqpos(cap_pose.residue(i), assembled.total_residue(),0)
        #print assembled.sequence()
        #print "LinkerB"
        for i in range(1, linkerB_pose.total_residue()+1):
            assembled.append_polymer_residue_after_seqpos(linkerB_pose.residue(i), assembled.total_residue(),0)
        #print assembled.sequence()
        #print "Bot2"
        for i in range( int(loodo_ins_begin)+3, bot.total_residue()+1): 
            assembled.append_polymer_residue_after_seqpos(bot.residue(i), assembled.total_residue(), 0)

        assembled.dump_pdb('Scaffold_{}.pdb'.format(base_name))

    

if __name__ == "__main__":
    main(sys.argv[1:])
