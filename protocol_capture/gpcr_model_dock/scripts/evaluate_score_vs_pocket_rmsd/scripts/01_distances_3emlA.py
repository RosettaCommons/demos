#################################################################################
# Calling:
# 	pymol -c 01_filter.py
#
# The script requires that the infilename (for example 1u19A.lst) is set  
# manually in the script. (the pymol interpreter is not fond of flags).
#
# The infilename contain a list of pdb file names (without .pdb)
# which are located in structure/infilename.split('.')[0]/ 
#
# It use the crystal structures in "crystal_pockets/" as templates for 
# calculating the distances.
#################################################################################

cmd.do('rei')
import sys, os
import glob
# from time import gmtime, strftime
# f_time = open('time.txt','w')
# f_time.write(strftime("%Y-%m-%d %H:%M:%S", gmtime())+'\n')

# Define list of crystal pockets and working dir.
infilename = '3emlA.lst' # 1u19A.lst

pdbid_to_be_filtered = infilename.split('.')[0]
working_dir = os.getcwd() + '/'
crystal_dir = working_dir + 'crystal_pockets/'
model_dir = working_dir + 'structures/'+pdbid_to_be_filtered+'/'
crystal_pocket_list_all =  ['1u19_rhodopsin_pocket','2rh1_b2ar_pocket','2vt4_b1ar_pocket','3eml_a2a_pocket','3odu_cxcr4_pocket','3pbl_dopamine_pocket','3rze_H1_pocket','3uon_M2_pocket','3v2w_sphingosine_pocket','4daj_m3_pocket','4djh_ka_opioid_pocket','4dkl_mu_opioid_pocket','4ea3_nociceptin_pocket','4ej4_de_opioid_pocket']
cmd.do("cd %s" % (working_dir))

# Avoid filtering based on the structure it self:
crystal_pocket_list = []
for structure in crystal_pocket_list_all:
	if not pdbid_to_be_filtered[0:4] in structure: crystal_pocket_list.append(structure)
print crystal_pocket_list

# Which residues are pocket residues in this pdb?
print 'here1 ' + crystal_dir
print 'here2 ' + pdbid_to_be_filtered[0:4]

crystal_pocket_path = ''.join(glob.glob(crystal_dir+pdbid_to_be_filtered[0:4]+'*'))
pocket_file_name = open(crystal_pocket_path,'r')
residues_in_pocket = []
for line in pocket_file_name.readlines():
	if "CA" in line.split(): residues_in_pocket.append(line.split()[5])

#################################################################################
# Prepare crystal pockets														#							
#################################################################################

def unique(items):
    found = set([])
    keep = []
    for item in items:
        if item not in found:
            found.add(item)
            keep.append(item)
    return keep

# Load pockets
for struc in crystal_pocket_list:
	print struc
	cmd.do("load %s%s.pdb" % (crystal_dir,struc))

# Align all structures to 2rh1
cmd.do('load %stemplate_2rh1.pdb' % (crystal_dir))
for pocket in crystal_pocket_list:
	cmd.do("cealign template_2rh1 and i. 82+86+90+93+106+109+110+113+114+117+118+121+199+200+203+204+208+282+286+289+290+292+293+296+309+311+312+315+316, %s" % (pocket))

# Delete all but ca-atoms
cmd.do('select name ca')
cmd.select("sele","(not sele)",enable=1)
cmd.remove("sele");cmd.delete("sele")

# Renumber all crystal pockets
for pocket in crystal_pocket_list:
	residue_list = []
	cmd.do("iterate %s,residue_list.append(resi)" % pocket) 
	residue_list_short = unique(residue_list)
	n = 1
	for residue_no in residue_list_short:		
		cmd.do("alter %s and resi %s,resi=%s" % (pocket,residue_no,n))
		n += 1
 
################################################################################
# Calculate all minimum distances between target pocket resides  			   #
# and all equivalent pocket residues in crystal structures       			   #							
################################################################################
 
# For which pdbs should distances be calculated?
with open('pdb_lists/'+infilename,'r') as f:
	pdb_list = f.read().splitlines()

distance_matrix = {}
for structure in pdb_list:
	distance_matrix[structure]=[]
	cmd.do("load %s%s.pdb" % (model_dir,structure))
	# Keep only pocket ca atoms
	cmd.do("remove %s and not resi %s" % (structure,'+'.join(residues_in_pocket)))
	cmd.do("remove %s and not name ca" % (structure))
	
	cmd.do("cealign template_2rh1 and i. 82+86+90+93+106+109+110+113+114+117+118+121+199+200+203+204+208+282+286+289+290+292+293+296+309+311+312+315+316, %s" % (structure))

	#renumber
	residue_list = []
	cmd.do("iterate %s,residue_list.append(resi)" % structure) 
	residue_list = unique(residue_list)
	n = 1
	for residue_no in residue_list:
		cmd.do("alter %s and resi %s,resi=%s" % (structure,residue_no,n))
		n += 1
	
	# calculate distances
	for residue_no in range(1,30):
		dists = []
		for crystal_pocket in crystal_pocket_list:
			if 'ca_' + crystal_pocket != structure: 
				dists.append(round(cmd.distance('dist','%s///%s/ca' % (structure, residue_no),'%s///%s/ca' % (crystal_pocket, residue_no)),3))
		distance_matrix[structure].append(min(dists))
		#print distance_matrix[structure]
		cmd.do("delete dist")
	cmd.do("delete %s" % structure)

print distance_matrix

#Filter all pdbs models
f = open('distance_output/'+pdbid_to_be_filtered+'_pymol_distance_output.txt','w')

for k in distance_matrix:
	distances = [k, ' '.join(map(str,distance_matrix[k]))]
  	for item in distances:
  		f.write(item+'\t')
  	f.write('\n')

# clean up
f.close()

# f_time.write(strftime("%Y-%m-%d %H:%M:%S", gmtime())+'\n')
