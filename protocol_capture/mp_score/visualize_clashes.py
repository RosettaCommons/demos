#! /usr/bin/python

# Script that creates a pymol input script to visualize scores
# fa_rep and dunbrack scores, but others can be set too 
# usage: 
#   visualize_clashes.py -p <input.pdb>
#   this writes a 'visualize_clashes.pml' script
#   open decoy in Pymol
#   run this script in pymol with command '@visualize_clashes.pml' 

import math, commands, os, sys, shutil
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def main(args):

  parser = OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
  parser.set_description(main.__doc__)

  parser.add_option('--pdb', '-p',
    action="store", default=None,
    help="PDB input file with MEM coordinates.", )

  parser.add_option('--rep', '-r',
    action="store_true", default=None,
    help="Use fa_repulsion to visualize. Default.", )

  parser.add_option('--dun', '-d',
    action="store_true", default=None,
    help="Use Dunbrack term to visualize.", )

  #parse options
  (options, args) = parser.parse_args(args=args[1:])
  global Options
  Options = options
	
	# output
  output = "visualize_clashes.pml"

  # Read in scores
  if Options.dun:
    score = (commands.getoutput("grep -E '^ALA|^CYS|^ASP|^GLU|^PHE|^GLY|^HIS|^ILE|^LYS|^LEU|^MET|^ASN|^PRO|^GLN|^ARG|^SER|^THR|^VAL|^TRP|^TYR' " + Options.pdb + " | awk '{print $18}'") ).splitlines()
    print "colored by Dunbrack score"
  else:
    score = (commands.getoutput("grep -E '^ALA|^CYS|^ASP|^GLU|^PHE|^GLY|^HIS|^ILE|^LYS|^LEU|^MET|^ASN|^PRO|^GLN|^ARG|^SER|^THR|^VAL|^TRP|^TYR' " + Options.pdb + " | awk '{print $3}'") ).splitlines()
#    score = (commands.getoutput("grep RSD_WTD " + Options.pdb + " | awk '{print $6}'") ).splitlines()
    print "colored by fa_rep"
	
  score = map( float, score )
  print score
		
  with file( output, 'w') as f:
    f.write("bg white\n")
    f.write("show cartoon\n")
    f.write("color silver\n")

    f.write("#ray shadows\n")
    f.write("set ray_shadows=0\n\n")

    # go through scores
    for i in range(0, len(score)):
			
#      color = "silver"
			
#      if score[i] <= 1.0:
#        color = "gray60"
      if score[i] > 1.0 and score[i] <= 10.0:
        color = "yellow"
        f.write("color " + color + ", resi " + str(i+1) + "\n")
      elif score[i] > 10.0 and score[i] <= 100.0:
        color = "orange"
        f.write("color " + color + ", resi " + str(i+1) + "\n")
      elif score[i] > 100.0:
        color = "red"
        f.write("color " + color + ", resi " + str(i+1) + "\n")

 
  print "wrote", output

################################################################################


if __name__ == "__main__" : main(sys.argv)

