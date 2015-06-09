The resfile contains reversions back to the native sequence of the
protein. The script generates a template residue file where one can
change residue back to native or test certain residues at different
position
	python generate_residuefile.py rosetta_cst.pdb
returns a resfile with type NATAA, NATRO (keep rotamer) - one should
add PIKAA for picking new residue e.g.
    1 PIKAA WFY 
