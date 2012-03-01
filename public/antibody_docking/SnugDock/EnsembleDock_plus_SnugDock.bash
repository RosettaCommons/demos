#!/bin/bash
/Users/jianqing/Work/Code_Development/rosetta-2.3.1/rosetta++/rosetta.mactel aa ABRM _ -dock -dock_mcm -quiet -nstruct 2 -fake_native -fab1 -pose -ensemble1 10 -dock_pert 3 8 8 -spin -ex1 -ex2aro_only -unboundrot -s ABRM -dock_rtmin -find_disulf -norepack_disulf -use_pdb_numbering -fake_native -skip_missing_residues -pose -snugdock -snugloop -snugh3 -snugh2
