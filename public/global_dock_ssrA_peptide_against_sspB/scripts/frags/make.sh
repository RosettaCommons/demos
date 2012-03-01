#!/usr/bin/tcsh

set PATH_TO_EXE=/vol/ek/ravehb/rosetta/svn_mini/rosetta_source/bin/
set PATH_TO_DB=/vol/ek/ravehb/rosetta/svn_mini/rosetta_database/

$PATH_TO_EXE/fragment_picker.linuxgccrelease -database $PATH_TO_DB @flags 

echo
