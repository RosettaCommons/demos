#!/usr/bin/tcsh

set PATH_TO_EXE=/vol/ek/liorz06/miniRosettaWorkspace/mini/bin/
set PATH_TO_DB=/vol/ek/liorz06/miniRosettaWorkspace/rosetta_database/

$PATH_TO_EXE/fragment_picker.linuxgccrelease -database $PATH_TO_DB @flags 

