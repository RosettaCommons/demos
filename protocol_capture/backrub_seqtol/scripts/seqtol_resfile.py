#!/usr/bin/python

import math
import os
import re
import sys
import types
import UserDict

residue_dict = {"ALA": "ALA", "CYS": "CYS", "ASP": "ASP", "ASH": "ASP", 
                "GLU": "GLU", "GLH": "GLU", "PHE": "PHE", "GLY": "GLY", 
                "HIS": "HIS", "HIE": "HIS", "HIP": "HIS", "ILE": "ILE", 
                "LYS": "LYS", "LYN": "LYS", "LEU": "LEU", "MET": "MET", 
                "ASN": "ASN", "PRO": "PRO", "GLN": "GLN", "ARG": "ARG", 
                "ARN": "ARG", "SER": "SER", "THR": "THR", "VAL": "VAL",
                "TRP": "TRP", "TYR": "TYR", "MSE": "MET"}

residues = ["ALA", "CYS", "ASP", "ASH", "GLU", "GLH", "PHE", "GLY", "HIS", 
            "HIE", "HIP", "ILE", "LYS", "LYN", "LEU", "MET", "ASN", "PRO", 
            "GLN", "ARG", "ARN", "SER", "THR", "VAL", "TRP", "TYR"]

canonical_residues = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", 
                      "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", 
                      "THR", "VAL", "TRP", "TYR"]

aa1_to_aa3 = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", 
              "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", 
              "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG", 
              "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}

aa3_to_aa1 = {"ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", 
              "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", 
              "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R", 
              "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y"}


class PDB:
	"""A class to store and manipulate PDB data"""
	
	def __init__(self, pdb = None):
		
		self.lines = []
		if type(pdb) == types.StringType:
			self.read(pdb)
		elif type(pdb) == types.ListType:
			self.lines.extend(pdb)
		elif type(pdb) == type(self):
			self.lines.extend(pdb.lines)
	
	def read(self, pdbpath):
		
		pdbhandle = open(pdbpath)
		self.lines = pdbhandle.readlines()
		pdbhandle.close()
	
	def write(self, pdbpath):
		
		pdbhandle = open(pdbpath, "w")
		pdbhandle.write("".join(self.lines))
		pdbhandle.close()
		
	def remove_nonbackbone_atoms(self, resid_list):
		
		backbone_atoms = set(["N", "CA", "C", "O", "OXT"])
		
		resid_set = set(resid_list)
	
		self.lines = [line for line in self.lines if line[0:4] != "ATOM" or 
		                                             line[21:27] not in resid_set or
		                                             line[12:16].strip() in backbone_atoms]

	def fix_backbone_occupancy(self):
	
		backbone_atoms = set(["N", "CA", "C", "O"])
		
		for i in xrange(len(self.lines)):
			line = self.lines[i]
			if line.startswith("ATOM") and line[12:16].strip() in backbone_atoms and float(line[54:60]) == 0:
				self.lines[i] = line[:54] + "  1.00" + line[60:]

	def remove_hetatm(self):
	
		self.lines = [line for line in self.lines if not line.startswith("HETATM")]

	def canonize_residue_names(self):
	
		for i in xrange(len(self.lines)):
			line = self.lines[i]
			if line.startswith("ATOM") or line.startswith("HETATM"):
				self.lines[i] = line[:17] + residue_dict[line[17:20]] + line[20:]

	def aa_resids(self):
		
		atomlines = [line for line in self.lines if line[0:4] == "ATOM" and 
		                                            line[17:20] in residues]
		
		resid_set = set()
		resid_list = []
		
		for line in atomlines:
			if line[21:27] not in resid_set:
				resid_set.add(line[21:27])
				resid_list.append(line[21:27])
		
		return resid_list
	
	def remove_resids(self, resids):
	
		resid_set = set(resids)
		self.lines = [line for line in self.lines if ((not (line.startswith("ATOM") or line.startswith("HETATM"))) or (line[21:27] not in resid_set))]
	
	def sort_chains(self):
	
		def chain_cmp(x, y):
			if (x.startswith("ATOM") or x.startswith("HETATM")) and (y.startswith("ATOM") or y.startswith("HETATM")):
				if x[21] < y[21]:
					return -1
				if x[21] > y[21]:
					return 1
			return 0
				
		self.lines.sort(cmp = chain_cmp)
	
	def atomlines(self, resid_list = None):
	
		if resid_list == None:
			resid_list = self.aa_resids()
		
		resid_set = set(resid_list)
		
		return [line for line in self.lines if line[0:4] == "ATOM" and 
		                                       line[21:27] in resid_set]

	def neighbors(self, distance, atom = " CA ", resid_list = None):
		
		if type(atom) != type(re.compile("")):
			lines = [line for line in self.atomlines(resid_list) if line[12:16] == atom]
		else:
			lines = [line for line in self.atomlines(resid_list) if atom.match(line[12:16])]
		
		shash = SpatialHash(distance)
		
		resid_pos = []
		
		for line in lines:
			pos = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
			shash.insert(pos, line[21:27])
		
		neighbor_dict = {}
		
		for line in lines:
			
			resid = line[21:27]
			pos = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
			
			neighbor_dict.setdefault(resid, [])
			neighbor_dict[resid] += [data[1] for data in shash.nearby(pos, distance)]
		
		return neighbor_dict

class ChainSequences(UserDict.DictMixin):
	"""A class for holding PDB chain sequences"""
	
	def __init__(self):
		self.seqs = {}
		self.chains = []
	
	def __getitem__(self, key):
		if type(key) == types.IntType:
			return self.seqs[self.chains[key]]
		else:
			return self.seqs[key]
	
	def __setitem__(self, key, value):
		if type(key) == types.IntType:
			self.seqs[self.chains[key]] = value
		else:
			self.seqs[key] = value
			if key not in self.chains:
				self.chains += [key]
	
	def __delitem__(self, key):
		if type(key) == types.IntType:
			del self.seqs[self.chains[key]]
			del self.chains[key]
		else:
			del self.seqs[key]
			self.chains.remove(key)
		
	def keys(self):
		return self.chains
	
	def parse_seqres(self, pdb):
		"""Parse the SEQRES entries into the object"""
		
		seqresre = re.compile("SEQRES")
		
		seqreslines = [line for line in pdb.lines if seqresre.match(line)]
		
		for line in seqreslines:
			chain = line[11]
			resnames = line[19:70].strip()
			self.setdefault(chain, [])
			self[chain] += resnames.split()

	def parse_atoms(self, pdb):
		"""Parse the ATOM entries into the object"""
		
		atomre = re.compile("ATOM")
		
		atomlines = [line for line in pdb.lines if atomre.match(line)]
		
		chainresnums = {}
		
		for line in atomlines:
			chain = line[21]
			resname = line[17:20]
			resnum = line[22:27]
			chainresnums.setdefault(chain, [])
			
			if resnum in chainresnums[chain]:
				assert self[chain][chainresnums[chain].index(resnum)] == resname
			else:
				self.setdefault(chain, [])
				self[chain] += [resname]
				chainresnums[chain] += [resnum]
		
		return chainresnums

	def seqres_lines(self):
		"""Generate SEQRES lines representing the contents"""
	
		lines = []
	
		for chain in self.keys():
			seq = self[chain]
			serNum = 1
			startidx = 0
			while startidx < len(seq):
				endidx = min(startidx+13, len(seq))
				lines += ["SEQRES  %2i %s %4i  %s\n" % (serNum, chain, len(seq), " ".join(seq[startidx:endidx]))]
				serNum += 1
				startidx += 13
		
		return lines
	
	def replace_seqres(self, pdb, update_atoms = True):
		"""Replace SEQRES lines with a new sequence, optionally removing 
		mutated sidechains"""
	
		newpdb = PDB()
		inserted_seqres = False
		entries_before_seqres = set(["HEADER", "OBSLTE", "TITLE",  "CAVEAT", "COMPND", "SOURCE", 
		                             "KEYWDS", "EXPDTA", "AUTHOR", "REVDAT", "SPRSDE", "JRNL", 
		                             "REMARK", "DBREF",  "SEQADV"])
	
		mutated_resids = {}
		
		if update_atoms:
			old_seqs = ChainSequences()
			chainresnums = old_seqs.parse_atoms(pdb)
		
			assert self.keys() == old_seqs.keys()
			
			for chain in self.keys():
				assert len(self[chain]) == len(old_seqs[chain])
				for i in xrange(len(self[chain])):
					if self[chain][i] != old_seqs[chain][i]:
						resid = chain + chainresnums[chain][i]
						mutated_resids[resid] = self[chain][i]
		
		for line in pdb.lines:
		
			entry = line[0:6]
			if (not inserted_seqres) and entry not in entries_before_seqres:
				inserted_seqres = True
				newpdb.lines += self.seqres_lines()
			
			if update_atoms and entry == "ATOM  ":
				resid = line[21:27]
				atom = line[12:16].strip()
				if not mutated_resids.has_key(resid):
					newpdb.lines += [line]
				else:
					newpdb.lines += [line[:17] + mutated_resids[resid] + line[20:]]
			elif entry != "SEQRES":
				newpdb.lines += [line]
		
		if update_atoms:
			newpdb.remove_nonbackbone_atoms(mutated_resids.keys())
		
		return newpdb
	
	def aa1_string(self, chain):
	
		return "".join([aa3_to_aa1[aa] for aa in self.seqs[chain]])

class SpatialHash:

	def __init__(self, size):
	
		self.size = size
		self.dimensions = 0
		self.quads = {}
	
	def quadkey(self, position):
	
		if len(position) != self.dimensions:
			sys.exit()
		
		quadrant = [0.]*self.dimensions
		for i in xrange(self.dimensions):
			quadrant[i] = int(math.floor(position[i]/self.size))
		
		return tuple(quadrant)
	
	def insert(self, position, data):
	
		if self.dimensions == 0:
			self.dimensions = len(position)

		key = self.quadkey(position)
		self.quads[key] = self.quads.get(key, []) + [(position, data)]

	def nearby(self, position, radius):
	
		minkey = self.quadkey([position[i] - radius for i in xrange(self.dimensions)])
		maxkey = self.quadkey([position[i] + radius for i in xrange(self.dimensions)])
		
		quadstosearch = [[i] for i in range(minkey[0], maxkey[0]+1)]
		
		for i in xrange(1, self.dimensions):
			newquads = []
			for j in xrange(minkey[i], maxkey[i]+1):
				newquads += [oldquad + [j] for oldquad in quadstosearch]
			quadstosearch = newquads
		
		quadstosearch = [tuple(quad) for quad in quadstosearch if self.quads.has_key(tuple(quad))]
		
		radiussquared = radius*radius
		
		results = []
		
		for quad in quadstosearch:
			for pos, data in self.quads[quad]:
				distsquared = 0
				for i in xrange(self.dimensions):
					distsquared += (position[i] - pos[i]) ** 2
				if distsquared <= radiussquared:
					results += [(pos, data)]
		
		return results

if __name__ == "__main__":

	if len(sys.argv) < 4:
	
		print "Usage: ./seqtol_resfile pdb_path design_command pos_1 [pos_2] [...]"
		sys.exit(1)

	pdb_path = sys.argv[1]
	
	pdb = PDB(pdb_path)

	resid_regexps = [res[0]+"[ ]*"+res[2:]+"[ ]*$" for res in sys.argv[2:]]
	
	resid_regexp = re.compile("("+"|".join(resid_regexps)+")")

	resids = pdb.aa_resids()
	
	resids_design = [resid for resid in resids if resid_regexp.match(resid)]
	
	distance = 10.
	atom = " CA "
	
	resids_neighbors = pdb.neighbors(distance, atom)
	resids_repack = set()
	
	for resid in resids_design:
		resids_repack.update(resids_neighbors[resid])
	
	seqtol_resfile_lines = ["NATRO\n", "start\n"]
	
	design_command = "ALLAA"
	design_command = "PIKAA ADEFGHIKLMNPQRSTVWY"
	design_command = sys.argv[2]
	
	for seqpos in xrange(len(resids)):
		resfile_id = resids[seqpos][1:].strip() + " " + resids[seqpos][0].strip()
		if resids[seqpos] in resids_design:
			seqtol_resfile_lines.append(resfile_id + " "+design_command+"\n")
		elif resids[seqpos] in resids_repack:
			seqtol_resfile_lines.append(resfile_id + " NATAA\n")
	
	base_path = os.path.splitext(pdb_path)[0]
	
	resfile_path = base_path+"_seqtol.resfile"
	
	seqtol_resfile = open(resfile_path, "w")
	seqtol_resfile.writelines(seqtol_resfile_lines)
	seqtol_resfile.close()
	
	seq = ChainSequences()
	seq.parse_atoms(pdb)
	
	seqs = []
	
	for chain in seq.chains:
		seqs += seq[chain]
	
	print "Output Resfile:"
	print resfile_path
	
	print "\nDesigned Residues:"
	for i in xrange(len(resids)):
		resid = resids[i]
		if resid in resids_design:
			print resid[0] + ":" + resid[1:].strip() + " " + seqs[i]
		
	print "\nRepacked Residues (" + str(distance) + " angstrom radius between " + atom.strip() + " atoms):"
	for i in xrange(len(resids)):
		resid = resids[i]
		if resid in resids_repack:
			print resid[0] + ":" + resid[1:].strip() + " " + seqs[i]
	
