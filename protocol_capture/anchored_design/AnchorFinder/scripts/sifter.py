#!/usr/bin/env python

#Args: first argument is integer cutoff; do not print lines below cutoff.

#The purpose of this script is to take "goodfile.out" style files from AnchorFinder, and sift out only those lines with cross-interface contacts with only one residue.

#For example: this line represents a single chain A to chain B interaction in a PDB with 4 chains.
#PDB pdb2vk1 window 112 loopness 5 nbrs 0 47 0 0 start 113 A pymol select pdb2vk1 and chain A and resi 113-117

#Assuming a cutoff of 40, t would pass the filter and be printed with the #neighbors duplicated to the front of the line for easy sorting:
#47 PDB pdb2vk1 window 112 loopness 5 nbrs 0 47 0 0 start 113 A pymol select pdb2vk1 and chain A and resi 113-117
#Assuming you pass a cutoff of 50, the example above would NOT be printed.

#This line has neighbors on multiple chains - it does not represent a dimeric interface.  It would not be printed.  (This is fake example data)
#PDB pdb2vk1 window 112 loopness 5 nbrs 0 47 0 54 start 113 A pymol select pdb2vk1 and chain A and resi 113-117


import string
import sys

if __name__ == '__main__':

    print "usage: sifter.py cutoff target_file"

    datafile = open(sys.argv[2])
        
    for eachline in datafile:

        just_nbrs = eachline[(eachline.find("nbrs ")+5):eachline.find(" start")]
        window_best = 0
        more_than_one = False
        for eachValue in just_nbrs.split():
            eachVal = int(eachValue)
            if ( not more_than_one and (window_best != 0) and (eachVal != 0) ) : more_than_one = True
            if (eachVal > window_best) : window_best = eachVal
        if ((window_best > int(sys.argv[1])) and not more_than_one):
            print window_best, eachline.rstrip()



