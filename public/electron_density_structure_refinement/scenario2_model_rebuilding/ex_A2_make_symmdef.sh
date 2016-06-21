#!/bin/bash

~/rosetta_workshop/rosetta/main/source/src/apps/public/symmetry/make_symmdef_file.pl \
 -m NCS -a A -i B C \
 -p setup_symm.pdb -r 1000 > D7.symm
