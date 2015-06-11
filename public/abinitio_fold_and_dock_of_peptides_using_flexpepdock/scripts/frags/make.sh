#!/usr/bin/csh

set PATH_TO_EXE=/vol/ek/ravehb/rosetta/svn_mini/mini/bin
set PATH_TO_DB=/vol/ek/ravehb/rosetta/svn_mini/minirosetta_database/

$PATH_TO_EXE/picker.linuxgccrelease -database $PATH_TO_DB @flags 

echo
