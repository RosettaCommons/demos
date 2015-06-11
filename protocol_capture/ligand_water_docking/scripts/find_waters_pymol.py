#!/sb/apps/Linux/bin/python

import pymol
pymol.finish_launching()
pymol.cmd.feedback('disable', 'all', 'actions')
pymol.cmd.feedback('disable', 'all', 'results')

#def add_interface_waters_to_pdb(complex_pdb, waters_pdb):
#	p=PDBParser(PERMISSIVE=1)
#	complex=p.get_structure('complex', complex_pdb)
#	waters=p.get_structure('complex', waters_pdb)
#	protein_residues=get_protein_residues(complex)
#	water_residues= waters.get_residues()
#	INH_atoms= complex[0]['X'][1].get_list()
#	water_residues_near_INH= filter(near_INH, water_residues)
#	water_residues_near_INH_and_protein= filter(near_protein, water_residues_near_INH)
#	
#	complex.add('W')
#	for water in water_residues_near_INH_and_protein:
#		complex['W'].add(water)
#
#	io=PDBIO()
#	io.set_structure(complex)
#	base=complex_pdb.partition('.')[0]
#	io.save( base+'with_interface_waters.pdb' )

from sys import exit
def has_interface_water_contacts(water_resi):
	distance_cutoff=2.5 # 3.0 for 'loose' and 'tight', 2.5 for 'extra_tight'
	inhibitor_contacts= pymol.cmd.find_pairs('waters AND resi '+water_resi, 'inhibitor', cutoff=distance_cutoff)
	protein_contacts= pymol.cmd.find_pairs('waters AND resi '+water_resi, 'protein', cutoff=distance_cutoff)
	#return len(inhibitor_contacts) > 0 and len(protein_contacts) > 0 # for 'loose_waters.pdb'
	return len(inhibitor_contacts) > 1 and len(protein_contacts) > 1 # for 'tight_waters.pdb'

def remove_non_interface_waters(protein, inhibitor, waters):
	pymol.cmd.load(protein, 'protein')
	pymol.cmd.load(inhibitor, 'inhibitor')
	pymol.cmd.load(waters, 'waters')
	myspace={'water_res_ids': []}
	pymol.cmd.iterate('waters AND resn WAT', 'water_res_ids.append(resi)', space=myspace)
	for water_resi in myspace['water_res_ids']:
		if not has_interface_water_contacts(water_resi):
			pymol.cmd.remove('waters and resi '+water_resi )
	pymol.cmd.save('extra_tight_waters.pdb', 'waters')
	pymol.cmd.quit()

from sys import argv
def main():
	args=argv[1:]
	if len(args) == 0:
		print 'Usage: find_waters_pymol.py protein.pdb inhibitor.pdb waters.pdb'
	assert(len(args)==3)
	protein= args[0]
	inhibitor= args[1]
	waters=args[2]
	remove_non_interface_waters(protein, inhibitor, waters)

if __name__ == '__main__':
	main()
