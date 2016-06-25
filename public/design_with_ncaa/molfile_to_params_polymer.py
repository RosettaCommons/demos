#!/usr/bin/env python
'''
Functions and executable for taking a ligand from an MDL Molfile
and splitting it into one or more .params files for Minirosetta.
See main() for usage or run with --help.

Author: Ian W. Davis
'''
import os, sys
if not hasattr(sys, "version_info") or sys.version_info < (2,4):
    raise ValueError("Script requires Python 2.4 or higher!")

from optparse import OptionParser
from sets import Set

# Magic spell to make sure Rosetta python libs are on the PYTHONPATH:
sys.path.append( os.path.dirname(os.path.dirname( os.path.abspath(sys.path[0]) )) )

from rosetta_py.io.mdl_molfile import *
from rosetta_py.utility.rankorder import argmin
from rosetta_py.utility import r3

# Features from Python 2.5 that we want to use:
if not hasattr(__builtins__, "any"):
    def any(itr):
        for el in itr:
            if el: return True
        return False

if not hasattr(__builtins__, "all"):
    def all(itr):
        for el in itr:
            if not el: return False
        return True


def add_fields_to_atoms(atoms):
    '''Adds a bunch of member variable placeholders that we use.'''
    for atom in atoms:
        atom.orig_name = atom.name # for kinemage output
        atom.pdb_name = ""         # PDB style atom name
        atom.ros_type = ""         # Rosetta atom type
        atom.mm_type = ""          # Molec. mechan. atom type
        atom.is_virtual = False    # for atom typing and charge assignment
        atom.rigid_id = 0          # non-zero id for regions with no rotatable bonds; may span fragments
        atom.fragment_id = 0       # non-zero id for fragments after bond breaking
        atom.conn_bonds = []       # list of cross-fragment bonds to this atom
        atom.is_root = False       # for atom tree
        atom.parent = None         # for atom tree
        atom.children = []         # for atom tree
        atom.stub1 = None          # for internal coords, derived from atom tree
        atom.stub2 = None          # for internal coords, derived from atom tree
        atom.stub3 = None          # for internal coords, derived from atom tree
        atom.input_stub1 = None    # for internal coords, derived from atom tree
        atom.input_stub2 = None    # for internal coords, derived from atom tree
        atom.input_stub3 = None    # for internal coords, derived from atom tree
        atom.d = 0.0               # distance to input_stub1
        atom.theta = 0.0           # 180 - angle with input_stub2 (degrees)
        atom.phi = 0.0             # dihedral from input_stub3 (degrees)
        atom.poly_upper = False    # is upper connect atom for polymer residue type
        atom.poly_lower = False    # is lower connect atom for polymer residue type
        atom.poly_n_bb = False     # is backbone nitrogen for polymer residue type
        atom.poly_ca_bb = False    # is backbone alpha carbon for polymer residue type
        atom.poly_c_bb = False     # is backbone carbonyl carbon for polymer residue type
        atom.poly_o_bb = False     # is backbone carbonyl oxygen for polymer residue type
        atom.poly_backbone = False # convience boolean
        atom.poly_ignore = False   # convience boolean
        atom.pdb_prefix_num = " "  # blah
        atom.pdb_elem = atom.elem  # blah
        atom.pdb_greek_dist = " "  # blah
        atom.pdb_postfix_num = " " # blah

def add_fields_to_bonds(bonds):
    '''Adds a bunch of member variable placeholders that we use.'''
    for bond in bonds:
        bond.can_rotate = False     # true for single bonds not in rings
        bond.is_proton_chi = False  # true for bonds that rotate OH hydrogens, etc
        bond.connection_id = 0      # non-zero id if bond crosses fragments
        # Remember we have to update mirror too!
        bond.mirror.can_rotate      = bond.can_rotate
        bond.mirror.is_proton_chi   = bond.is_proton_chi
        bond.mirror.connection_id   = bond.connection_id
        bond.poly_ignore = False    # convience boolean

def find_virtual_atoms(atoms):
    '''Atoms whose names start with "V" are virtual, used in enzyme design, etc.'''
    for atom in atoms:
        if atom.name.startswith("V") or atom.name.startswith("X"):
            atom.is_virtual = True
            atom.elem = "X" # "V" is a real element (Vanadium)

def uniquify_atom_names(atoms):
    '''If atom names are not unique, rename/number them.'''
    # First try to pad names to match PDB convention:
    for atom in atoms:
        # If two chars, assume it's an element name
        #if len(atom.name) >= 2 and atom.name[1].isalpha(): atom.name = "%-4s" % atom.name
        #else: atom.name = " %-3s" % atom.name
        # Assume that we got the element name right, and that the atom name starts with the element symbol
        if len(atom.elem) == 1 and len(atom.name) <= 3: atom.name = " %-3s" % atom.name
        else: atom.name = "%-4s" % atom.name
    duplicate_names = False
    atom_names = Set()
    for atom in atoms:
        if atom.name in atom_names:
            duplicate_names = True
            break
        atom_names.add(atom.name)
    if not duplicate_names: return
    print "Atom names contain duplications -- renaming all atoms."
    # This is potentially bad (> 4 char names) for > 100 atoms:
    #for i,atom in enumerate(atoms):
    #    atom.name = "%s%i" % (atom.elem, i+1)
    # So instead, number each element separately:
    atom_names.clear()
    for atom in atoms:
        i = 1
        while True:
            name = "%2s%i" % (atom.elem, i)
            if name not in atom_names: break
            i += 1
        atom.name = name
        atom_names.add(name)

def check_bond_count(atoms):
    '''Safety check for structures with stupid numbers of bonds to atoms, at Florian's request.'''
    def valence(bond):
        if bond.order == Bond.AROMATIC: return 1.5
        else: return float(bond.order)
    for atom in atoms:
        if atom.is_H and len(atom.bonds) > 1:
            raise ValueError("Atom %s is a hydrogen with >1 bonds" % atom.name)
        # Valence 4.5 for e.g. carbon comes up in peptide bonds and at the joins in multi-ring systems.
        if sum([valence(bond) for bond in atom.bonds]) > 4.5:
            print "WARNING:  atom %s has valence > 4" % atom.name

def check_aromaticity(bonds):
    '''Safety check for Kekule structures (alternating single/double bonds)
    rather than bonds described explicitly as aromatic.'''
    aro_bonds = len([b for b in bonds if b.order == Bond.AROMATIC])
    dbl_bonds = len([b for b in bonds if b.order == Bond.DOUBLE])
    if aro_bonds == 0 and dbl_bonds > 0:
        print "WARNING:  structure contains double bonds but no aromatic bonds"
        print "  Aromatic bonds must be identified explicitly --"
        print "  alternating single/double bonds (Kekule structure) won't cut it."
        print "  This warning does not apply to you if your molecule really isn't aromatic."

def assign_rosetta_types(atoms):
    '''Assigns Rosetta atom types.
    Based on Rosetta++ ligand_ns.cc set_rosetta_atom_types().
    This has been tested against the assignments produced by (Jens? Molecule.exe?)
    for the Meiler and Baker 2006 cross docking test set;
    cases where they disagree are usually due to weird
    bond orders (probably mistakes) in the .mol files.

    As I look through these rules, I see that they contain an ad-hoc attempt
    at aromaticity perception.  Not sure how general this is, though.
    '''
    # Helper function, count bonded atoms that match predicate
    # Predicate takes one arg, an atom
    def count_bonded(atom, pred):
        return sum(1 for bond in atom.bonds if not bond.a2.is_virtual and pred(bond.a2))
    # A predicate for count_bonded()
    def is_aromatic(atom):
        return any( bond.order == Bond.AROMATIC for bond in atom.bonds if not bond.a2.is_virtual )
    # For carbon and nitrogen, is_saturated implies SP3 hybridization
    # Not quite true:  look at ring N in 1aq1 -- 3 "single" bonds but flat (like Trp)
    def is_saturated(atom):
        return all( bond.order == Bond.SINGLE for bond in atom.bonds if not bond.a2.is_virtual )
    # For each atom, analyze bonding pattern to determine type
    for i, a in enumerate(atoms): # i just used for debugging output
        # H, C, O, N have complicated rules.
        # Everything else maps to a single atom type.
        if a.is_virtual:
            a.ros_type = "VIRT"
        elif "H"  == a.elem:
            num_aro_C = count_bonded(a, lambda x: (x.elem == "C" and is_aromatic(x)) or x.ros_type == "aroC")
            num_NOS = count_bonded(a, lambda x: x.elem == "N" or x.elem == "O" or x.elem == "S")
            if num_NOS >= 1:    a.ros_type = "Hpol"
            elif num_aro_C >=1: a.ros_type = "Haro"
            else:               a.ros_type = "Hapo"
        elif "C"  == a.elem:
            if is_saturated(a):
                num_H = count_bonded(a, lambda x: x.is_H)
                if num_H >= 3:      a.ros_type = "CH3 "
                elif num_H == 2:    a.ros_type = "CH2 "
                else:               a.ros_type = "CH1 "
            else:
                num_dbl_nonO = 0; num_aro_nonO = 0; num_aro_N = 0;
                a_bonds = [b for b in a.bonds if not b.a2.is_virtual]
                for bond in a_bonds:
                    if bond.order == Bond.DOUBLE:
                        if bond.a2.elem != "O": num_dbl_nonO += 1
                    elif bond.order == Bond.AROMATIC:
                        if bond.a2.elem != "O": num_aro_nonO += 1
                        if bond.a2.elem == "N": num_aro_N += 1 # really if, not elif
                #print i+1, a.name, num_aro_nonO, num_dbl_nonO, num_aro_N
                if num_aro_nonO >= 2:   a.ros_type = "aroC"
                elif num_dbl_nonO >= 1: a.ros_type = "aroC"
                elif num_aro_N >= 1:    a.ros_type = "CNH2"
                else:                   a.ros_type = "COO "
        elif "N"  == a.elem:
            num_H = count_bonded(a, lambda x: x.is_H)
            heavy_nbrs = count_bonded(a, lambda x: not x.is_H)
            assert num_H + heavy_nbrs == len(a.bonds)
            if num_H >= 3:              a.ros_type = "Nlys" # carries a VERY high desolvation penalty
            # Not totally sure about this one, may want Ntrp instead if more than one heavy neighbor:
            elif num_H == 2:            a.ros_type = "NH2O" # Narg would also be a possibility, but they're fairly similar
            elif num_H == 1:
                if heavy_nbrs <= 2:     a.ros_type = "Ntrp" # should always be 2 neighbors, not less
                else:                   a.ros_type = "Ntrp" # Npro? protonated tertiary amine
            else: # num_H == 0
                if heavy_nbrs <= 2:     a.ros_type = "Nhis"
                elif heavy_nbrs == 3:
                    if is_saturated(a): a.ros_type = "Nhis" # deprotonated tertiary amine; need an sp3 hybrid H-bond acceptor type...
                    # This also catches nitro groups -- is that what we want here?
                    else:               a.ros_type = "Npro" # X=[N+](X)X, including nitro groups
                else:                   a.ros_type = "Npro" # quaternary amine
        elif "O"  == a.elem:
            num_H = count_bonded(a, lambda x: x.is_H)
            num_bonds = count_bonded(a, lambda x: True) #len(a.bonds)  using count_bonded() avoid counting virtual bonds/atoms
            bonded_to_N = count_bonded(a, lambda x: x.elem == "N")
            bonded_to_C_to_N = count_bonded(a, lambda x: x.elem == "C" and count_bonded(x, lambda y: y.elem == "N") > 0)
            if is_saturated(a):
                if num_bonds >= 2:
                    unsat_nbrs = count_bonded(a, lambda x: not is_saturated(x))
                    if num_H > 0:   a.ros_type = "OH  " # catches C(=O)OH (Kekule form)
                    elif a.is_ring and a.ring_size < 5:
                                    a.ros_type = "OH  " # small, strained rings leave the O more exposed? (IWD, see 1p8d)
                    elif a.is_ring and unsat_nbrs > 0:
                                    a.ros_type = "Oaro" # catches aromatic O in furan-like rings, though I rarely see these H-bond (IWD)
                    else:           a.ros_type = "OH  " # catches ethers, ROR (IWD, see comment)
                    # The lone pairs on ethers are capable of H-bonding in the same way that alcohols are.
                    # While alkyl ethers are quite non-polar, many others seem to make Hbonds,
                    # such as those attached to phosphates (R-O-PO3), methyls (R-O-CH3), and aromatic rings (R-O-Ph).
                    # It is unclear from the literature how strong these are, and is probably very situation dependent.
                else:               a.ros_type = "OOC " # catches C(=O)[O-] (Kekule form) -- new rule by IWD
            elif num_H > 0:         a.ros_type = "OH  " # catches c(o)oH (aromatic bonds to both O)
            elif bonded_to_N:       a.ros_type = "ONH2"
            # This is a non-standard rule introduced by IWD, agreed to by KWK:
            elif bonded_to_C_to_N:  a.ros_type = "ONH2"
            else:                   a.ros_type = "OOC "
        elif "S"  == a.elem: a.ros_type = "S   "
        elif "P"  == a.elem: a.ros_type = "Phos"
        elif "F"  == a.elem: a.ros_type = "F   "
        elif "CL" == a.elem: a.ros_type = "Cl  "
        elif "BR" == a.elem: a.ros_type = "Br  "
        elif "I"  == a.elem: a.ros_type = "I   "
        elif "NA" == a.elem: a.ros_type = "Na1p"
        elif "K"  == a.elem: a.ros_type = "K1p "
        elif "MG" == a.elem: a.ros_type = "Mg2p"
        elif "FE" == a.elem: a.ros_type = "Fe3p"
        elif "CA" == a.elem: a.ros_type = "Ca2p"
        elif "ZN" == a.elem: a.ros_type = "Zn2p"
        else: raise ValueError("Unknown element '%s'" % a.elem)

def assign_mm_types(atoms, peptoid):
    '''Written by Doug Renfrew strongly influenced by the above function. This _TRYS_ to fill in the CHARMM27 atom types.
    It is rather conservative. It may but X in places in which case you will be tasked with finding the correct type, sorry.
    Rather than the crazy if/else tree (probably more appropriately forest) above I tried to split each type up in to
    its own function. There is unfortunatly some order dependence in calling all the functions but I have tried to avoid it. '''
    # helper functions
    def count_bonded(atom, pred):
        return sum(1 for bond in atom.bonds if not bond.a2.is_virtual and pred(bond.a2))
    def is_aromatic(atom):
        return any( bond.order == Bond.AROMATIC for bond in atom.bonds if not bond.a2.is_virtual )
    def is_saturated(atom):
        return all( bond.order == Bond.SINGLE for bond in atom.bonds if not bond.a2.is_virtual )
    # ----- HYDROGEN -----
    # ----- HYDROGEN -----
    # ----- HYDROGEN -----
    def is_charmm_H(atom, attached):
        ''' For polar hydogen: attached to O or one of 2 H attached to a N '''
        if attached.elem == "O": return True
        elif attached.elem == "N":
            attached_num_H =  count_bonded(attached, lambda x: x.is_H)
            if attached_num_H <= 2: return True
        else: return False
    def is_charmm_HC(atom, attached):
        ''' For nterm hydogen: one of 3 H attached to a N '''
        if attached.elem == "N":
            attached_num_H = count_bonded(attached, lambda x: x.is_H)
            if attached_num_H == 3: return True
        else: return False
    def is_charmm_HA(atom, attached):
        ''' For nonpolar hydrogen: an H attached to a saturated carbon '''
        attached_sat = count_bonded(atom, lambda x: (x.elem == "C" and is_saturated(x)))
        if attached_sat: return True
        else: return False
    def is_charmm_HP(atom, attached):
        ''' For aromatic hydrogen: hydrogen attached to an aromatic carbon '''
        num_aro_C = count_bonded(atom, lambda x: (x.elem == "C" and is_aromatic(x)))
        if num_aro_C >= 1: return True
    def is_charmm_HB(atom, attached, peptoid):
        ''' For backbone hydrogen: hydogen attached to a backbone carbon alpha '''
        if attached.poly_ca_bb and not peptoid: return True
        else: return False
    def is_charmm_HR1(atom, attached):
        ''' For nutral HIS HE1 hydrogen: (A) attached to an aromatic carbon that is between two N and only 
        one is prot or (B) attached to and aromatic carbon between an N and C and both ring N are protonated '''
        if is_aromatic(attached) and attached.elem == "C":
            attached_N = [bond.a2 for bond in attached.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "N"]
            attached_C = [bond.a2 for bond in attached.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "C"]
            # case A
            if len(attached_N) == 2:
                num_prot_N = sum(1 for an in attached_N if len(an.bonds) == 3 ) # catches H and CH3
                if num_prot_N == 1: return True
                else: return False
            # case B
            elif len(attached_N) == 1 and len(attached_C) == 1:
                attached_C_N = [bond.a2 for bond in attached_C[0].bonds if bond.order == Bond.AROMATIC and  bond.a2.elem == "N"]
                num_sat_N =  sum(1 for an  in attached_N   if len(an.bonds ) == 3 ) # catches H and CH3
                num_sat_N += sum(1 for acn in attached_C_N if len(acn.bonds) == 3 ) # catches H and CH3
                if num_sat_N == 2: return True
                else: return False
            else: return False
        else:return False
    def is_charmm_HR2(atom, attached):
        ''' For protonated HIS HE1: attached to aromatic C b/w two protonated N '''
        if is_aromatic(attached) and attached.elem == "C":
            attached_N = [bond.a2 for bond in attached.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "N"]
            if len(attached_N) == 2:
                num_prot_N = sum(1 for an in attached_N if len(an.bonds) == 3 ) # catches H and CH3
                if num_prot_N == 2: return True
                else: return False
            else: return False
        else:return False
    def is_charmm_HR3(atom, attached):
        ''' For nutral HIS HD2 hydrogen: attached to aromatic C b/w an N and C and one N on ring is protonated'''
        if is_aromatic(attached) and attached.elem == "C":
            attached_N = [bond.a2 for bond in attached.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "N"]
            attached_C = [bond.a2 for bond in attached.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "C"]
            if len(attached_N) == 2: return False
            elif len(attached_N) == 1 and len(attached_C) == 1:
                attached_C_N = [bond.a2 for bond in attached_C[0].bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "N"]
                num_sat_N =  sum(1 for an  in attached_N   if len(an.bonds ) == 3 ) # catches H and CH3
                num_sat_N += sum(1 for acn in attached_C_N if len(acn.bonds) == 3 ) # catches H and CH3
                if num_sat_N == 1: return True
                else: return False
            else: return False
        else:return False
    def is_charmm_HS(atom, attached):
        ''' For thiol hydrogen: hydrogen attached to a sulfer '''
        if attached.elem == "S": return True
    def is_charmm_HE1(atom, attached):
        ''' For alkene hydrogen: hydrogen of the type RHC=CR '''
        attached_sat = count_bonded(atom, lambda x: (x.elem == "C" and is_saturated(x)))
        attached_num_H = count_bonded(attached, lambda x: x.is_H)
        if not attached_sat and attached_num_H == 1: return True
        else: return False
    def is_charmm_HE2(atom, attached):
        ''' For alkene hydrogen: hydrogen of the type H2C=CR '''
        attached_sat = count_bonded(atom, lambda x: (x.elem == "C" and is_saturated(x)))
        attached_num_H = count_bonded(attached, lambda x: x.is_H)
        if not attached_sat and attached_num_H == 2: return True
        else: return False
    def is_charmm_HF1(atom, attached):
        ''' For Aliphatic H on fluorinated C: hydrogen attached to a carbon that is also attached to one fluorine'''
        if attached.elem == "C":
            attached_num_F = count_bonded(attached, lambda x: x.elem == "F")
            if attached_num_F == 1: return True
            else: return False
        else:return False
    def is_charmm_HF2(atom, attached):
        ''' For Aliphatic H on fluorinated C: hydrogen attached to a carbon that is also attached to two fluorine'''
        if attached.elem == "C":
            attached_num_F = count_bonded(attached, lambda x: x.elem == "F")
            if attached_num_F == 2: return True
            else: return False
        else:return False
    # ----- CARBON -----
    # ----- CARBON -----
    # ----- CARBON -----
    def is_charmm_CA(atom):
        ''' For aromatic carbon: carbon with 2 aromatic bonds '''
        num_hvy_bonds = sum(1 for bond in atom.bonds if not bond.a2.elem == "H")
        num_aro_bonds = sum(1 for bond in atom.bonds if bond.order == Bond.AROMATIC)
        if num_aro_bonds == 2 and num_hvy_bonds >= 2: return True
        else: return False
    def is_charmm_CT(atom):
        ''' For aliphatic sp3 C without hydrogens: carbon with all single bonds none to hydrogen '''
        if is_saturated(atom) and len(atom.bonds) == 4:
            num_H = count_bonded(atom, lambda x: x.elem == "H")
            if num_H == 0: return True
            else: return False
        else: return False
    def is_charmm_CT1(atom): 
        ''' For aliphatic sp3 C with 1 hydrogens: carbon with all single bonds one to hydrogen '''
        if is_saturated(atom) and len(atom.bonds) == 4:
            num_H = count_bonded(atom, lambda x: x.elem == "H")
            if num_H == 1: return True
            else: return False
        else: return False
    def is_charmm_CT2(atom):
        ''' For aliphatic sp3 C with 2 hydrogens: carbon with all single bonds two to hydrogen '''
        if is_saturated(atom) and len(atom.bonds) == 4:
            num_H = count_bonded(atom, lambda x: x.elem == "H")
            if num_H == 2: return True
            else: return False
        else: return False
    def is_charmm_CT3(atom):
        ''' For aliphatic sp3 C with 3 hydrogens: carbon with all single bonds three to hydrogen '''
        if is_saturated(atom) and len(atom.bonds) == 4:
            num_H = count_bonded(atom, lambda x: x.elem == "H")
            if num_H == 3: return True
            else: return False
        else: return False
    def is_charmm_CPT(atom):
        ''' For bridging carbons: carbon with 3 aromatic and 3 bonds to heavy atoms that bridges more than one aromatic ring'''
        num_hvy_bonds = sum(1 for bond in atom.bonds if not bond.a2.elem == "H")
        num_aro_bonds = sum(1 for bond in atom.bonds if bond.order == Bond.AROMATIC)
        if num_aro_bonds == 3 and num_hvy_bonds == 3: return True
        else: return False
    def is_charmm_CS(atom):
        ''' For thiolate carbon: carbon attached to a "bare" sulfer ie. H3C*S(-) '''
        attached_S = [bond.a2 for bond in atom.bonds if bond.a2.elem == "S"]
        if len(attached_S) == 1:
            if len(attached_S[0].bonds) == 1: return True # does the sulfer only make one bond
            else: return False
        else: return False
    def is_charmm_CE1(atom):
        ''' For alkene carbon RHC*=CR: carbon double bonded to another carbond and single bonded to a hydorgen and something else '''
        attached_double = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C" and bond.order == Bond.DOUBLE]
        num_H = count_bonded(atom, lambda x: x.is_H)
        if not is_saturated(atom) and len(attached_double) == 1 and num_H <= 1: return True
        else: return False
    def is_charmm_CE2(atom):
        ''' For alkene carbon H2C*=CR: carbon double bonded to another carbond and single bonded to two hydorgens '''
        attached_double = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C" and bond.order == Bond.DOUBLE]
        num_H = count_bonded(atom, lambda x: x.is_H)
        if not is_saturated(atom) and len(attached_double) == 1 and num_H == 2: return True
    def is_charmm_CN(atom):
        ''' For carbon in cyano (nitrial) group: carbon triple bonded to a nitrogen '''
        num_trip_N = sum(1 for bond in atom.bonds if bond.order == Bond.TRIPLE and bond.a2.elem == "N")
        if num_trip_N == 1: return True
        else: return False
    def is_charmm_CPH1(atom):
        ''' For imidazole carbon: carbon in imidazol bonded to a N and a C ie. his CG and CD2 carbons '''
        if is_aromatic(atom):
            attached_N = [bond.a2 for bond in atom.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "N"]
            attached_C = [bond.a2 for bond in atom.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "C"]
            if len(attached_N) == 1 and len(attached_C) == 1:
                attached_C_N = [bond.a2 for bond in attached_C[0].bonds if bond.a2.elem == "N"]
                if len(attached_C_N) == 1: return True
                else: return False
            else: return False
        else: return False
    def is_charmm_CPH2(atom):
        ''' For imidazole carbon: carbon in imidiazole bonded to two N ie. his CE1 carbon '''
        if is_aromatic(atom):
            attached_N = [bond.a2 for bond in atom.bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "N"]
            if len(attached_N) == 2: return True
            else: return False
        else: return False
    def is_charmm_CY(atom):
        ''' For indol carbon: carbon with aro bond to a bridge carbon (atom type CPT) '''
        attached_C = [bond.a2 for bond in atom.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "C"]
        for ac in attached_C:
            num_hvy_bonds = sum(1 for bond in ac.bonds if not bond.a2.elem == "H")
            num_aro_bonds = sum(1 for bond in ac.bonds if bond.order == Bond.AROMATIC)
            if num_aro_bonds == 3 and num_hvy_bonds == 3 and atom.ring_size == 5: return True
        return False
    def is_charmm_CC(atom):
        ''' For carboyl carbon: carbon double bonded to O and single bonded to an unprotonated O.
        Catches esters or carboxylic acids but not aldehyde, ketone, amides '''
        attached_double = [bond.a2 for bond in atom.bonds if bond.a2.elem == "O" and bond.order == Bond.DOUBLE]
        attached_single = [bond.a2 for bond in atom.bonds if bond.a2.elem == "O" and bond.order == Bond.SINGLE ]
        if len(attached_double) == 1 and len(attached_single) == 1:
            if len(attached_single[0].bonds) == 1: return True # O is only one bond to the carbon
            else: return False
        else:return False
    def is_charmm_CD(atom):
        ''' For carboyl carbon: carbon double bonded to O and single bonded to a protonated O.
        Catches esters or carboxylic acids but not aldehyde, ketone, amides '''
        attached_double = [bond.a2 for bond in atom.bonds if bond.a2.elem == "O" and bond.order == Bond.DOUBLE]
        attached_single = [bond.a2 for bond in atom.bonds if bond.a2.elem == "O" and bond.order == Bond.SINGLE ]
        if len(attached_double) == 1 and len(attached_single) == 1:
            if len(attached_single[0].bonds) > 1: return True # O is bonded to something in addition to the carbon
            else: return False
        else:return False
    def is_charmm_C(atom):
        ''' For carboyl carbon: carbon double bonded to N or O and single bonded to two other atoms that are not O.
        Catches aldehyde, ketone, amides but not esters or carboxylic acids which are a special atom type '''
        attached_double = [bond.a2 for bond in atom.bonds if (bond.a2.elem == "O" or bond.a2.elem == "N") and bond.order == Bond.DOUBLE]
        attached_single = [bond.a2 for bond in atom.bonds if not bond.a2.elem == "O" and bond.order == Bond.SINGLE ]
        if len(attached_double) == 1 and len(attached_single) == 2: return True
        else:return False
    # ----- NITROGEN -----
    # ----- NITROGEN -----
    # ----- NITROGEN -----
    def is_charmm_N(atom):
        ''' For proline nitrogen: nitrogen with 3 heavy bonds '''
        num_hvy_bond = sum(1 for bond in atom.bonds if not bond.a2.elem == "H")
        num_H = sum(1 for bond in atom.bonds if bond.a2.elem == "H")
        if num_hvy_bond == 3 and num_H == 0: return True
        else: return False
    def is_charmm_NR1(atom):
        ''' For neutral his protonated ring nitrogen: protonated aromatic nitrogen '''
        if is_aromatic(atom) and atom.ring_size == 5:
            attached_aro_C = [bond.a2 for bond in atom.bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "C"]
            attached_single = [bond.a2 for bond in atom.bonds if bond.order == Bond.SINGLE] # catches H and CH3
            if len(attached_aro_C) == 2 and len(attached_single) == 1:
                attached_aro_C_aro_N  = [bond.a2 for bond in attached_aro_C[0].bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "N" and bond.a2.name != atom.name]
                attached_aro_C_aro_N += [bond.a2 for bond in attached_aro_C[1].bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "N" and bond.a2.name != atom.name]
                if len(attached_aro_C_aro_N) == 1:
                    single = [bond.a2 for bond in attached_aro_C_aro_N[0].bonds if bond.order == Bond.SINGLE] # catches H and CH3
                    if len(single) == 0: return True
                    else: return False
                else: return False
        else:return False
    def is_charmm_NR2(atom):
        ''' For neutral his unprotonated ring nitrogen: unprotonated aromatic nitrogen '''
        if is_aromatic(atom) and atom.ring_size == 5:
            attached_aro_C = [bond.a2 for bond in atom.bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "C"]
            attached_single = [bond.a2 for bond in atom.bonds if bond.order == Bond.SINGLE] # catches H and CH3
            if len(attached_aro_C) == 2 and len(attached_single) == 0:
                attached_aro_C_aro_N  = [bond.a2 for bond in attached_aro_C[0].bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "N" and bond.a2.name != atom.name]
                attached_aro_C_aro_N += [bond.a2 for bond in attached_aro_C[1].bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "N" and bond.a2.name != atom.name]
                if len(attached_aro_C_aro_N) == 1:
                    single = [bond.a2 for bond in attached_aro_C_aro_N[0].bonds if bond.order == Bond.SINGLE] # catches H and CH3
                    if len(single) == 1: return True
                    else: return False
                else: return False
        else:return False
    def is_charmm_NR3(atom):
        ''' For charged his ring nitrogen: protonated nitrogen attached to a carbon attachec to another protonated nitrogen'''
        if is_aromatic(atom) and atom.ring_size == 5:
            attached_aro_C = [bond.a2 for bond in atom.bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "C"]
            attached_single = [bond.a2 for bond in atom.bonds if bond.order == Bond.SINGLE] # catches H and CH3
            if len(attached_aro_C) == 2 and len(attached_single) == 1:
                attached_aro_C_aro_N  = [bond.a2 for bond in attached_aro_C[0].bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "N" and bond.a2.name != atom.name]
                attached_aro_C_aro_N += [bond.a2 for bond in attached_aro_C[1].bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "N" and bond.a2.name != atom.name]
                if len(attached_aro_C_aro_N) == 1:
                    single = [bond.a2 for bond in attached_aro_C_aro_N[0].bonds if bond.order == Bond.SINGLE] # catches H and CH3
                    if len(single) == 1: return True
                    else: return False
                else: return False
        else:return False
    def is_charmm_NH1(atom):
        ''' For peptide nitrogen: nitrogen with 2 heavy bonds and 1 hydogen'''
        num_hvy_bond = sum(1 for bond in atom.bonds if not bond.a2.elem == "H")
        num_H = sum(1 for bond in atom.bonds if bond.a2.elem == "H")
        if num_hvy_bond == 2 and num_H == 1: return True
        else: return False
    def is_charmm_NH2(atom):
        ''' For amide nitrogen: nitrogen with 1 heavy bond and 2 hydogen'''
        num_hvy_bond = sum(1 for bond in atom.bonds if not bond.a2.elem == "H")
        num_H = sum(1 for bond in atom.bonds if bond.a2.elem == "H")
        if num_hvy_bond == 1 and num_H == 2: return True
        else: return False
    def is_charmm_NH3(atom):
        ''' For ammonium nitrogen: nitrogen with 1 heavy bond and 3 hydogen'''
        num_hvy_bond = sum(1 for bond in atom.bonds if not bond.a2.elem == "H")
        num_H = sum(1 for bond in atom.bonds if bond.a2.elem == "H")
        if num_hvy_bond == 1 and num_H == 3: return True
        else: return False
    def is_charmm_NC2(atom):
        ''' For guanidinium nitrogen: nitrogen attached to carbon that is bonded to two other nitrogens'''
        attached_C = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C"]
        if len(attached_C) == 1:
            attached_C_N = [bond.a2 for bond in attached_C[0].bonds if bond.a2.elem == "N"]
            if len(attached_C_N) == 3: return True
            else:return False
        else: return False
    def is_charmm_NY(atom):
        ''' For indol nitrogen: nitrogen with aro bond to a bridge carbon (atom type CPT) '''
        attached_C = [bond.a2 for bond in atom.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "C"]
        for ac in attached_C:
            num_hvy_bonds = sum(1 for bond in ac.bonds if not bond.a2.elem == "H")
            num_aro_bonds = sum(1 for bond in ac.bonds if bond.order == Bond.AROMATIC)
            if num_aro_bonds == 3 and num_hvy_bonds == 3 and atom.ring_size == 5: return True
        return False
    def is_charmm_NP(atom):
        ''' For nterm proline nitrogen: nitrogen with 2 heavy bond and 2 hydogen'''
        num_hvy_bond = sum(1 for bond in atom.bonds if not bond.a2.elem == "H")
        num_H = sum(1 for bond in atom.bonds if bond.a2.elem == "H")
        if num_hvy_bond == 2 and num_H == 2: return True
        else: return False
    def is_charmm_NC(atom):
        ''' For carbon in cyano (nitrial) group: carbon triple bonded to a nitrogen '''
        num_trip_C = sum(1 for bond in atom.bonds if bond.order == Bond.TRIPLE and bond.a2.elem == "C")
        if num_trip_C == 1: return True
        else: return False
    # ----- OXYGEN -----
    # ----- OXYGEN -----
    # ----- OXYGEN -----
    def is_charmm_O(atom):
        ''' For carbonyl oxygen: oxygen only double bonded to a carbon that is not bonded to another oxygen '''
        attached_double_C = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C" and bond.order == Bond.DOUBLE]
        attached_single_C = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C" and bond.order == Bond.SINGLE]
        if len(attached_double_C) == 1 and len(attached_single_C) == 0:
            attached_double_C_O = [bond.a2 for bond in attached_double_C[0].bonds if bond.a2.elem == "O" and bond.a2.name != atom.name]
            if len(attached_double_C_O) == 0: return True
            else: return False
        else: return False
    def is_charmm_OB(atom):
        ''' For carbonyl oxygen in acetic acid: oxygen only double bonded to a carbon that is single bonded to another oxygen that itself has a bond to something else'''
        attached_double_C = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C" and bond.order == Bond.DOUBLE]
        attached_single_C = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C" and bond.order == Bond.SINGLE]
        if len(attached_double_C) == 1 and len(attached_single_C) == 0:
            attached_double_C_O = [bond.a2 for bond in attached_double_C[0].bonds if bond.a2.elem == "O" and bond.a2.name != atom.name]
            if len(attached_double_C_O) == 1:
                attached_double_C_O_any = [bond.a2 for bond in attached_double_C_O[0].bonds if bond.a2.name != attached_double_C[0].name]
                if len(attached_double_C_O_any) > 0: return True
                else: return False
            else: return False
        else: return False
    def is_charmm_OC(atom):
        ''' For carboxylate oxygen: oxygen only (A)single/(B)double bonded to a carbon that is (A)double/(B)single bonded to another bare oxygen that is '''
        attached_double_C = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C" and bond.order == Bond.DOUBLE]
        attached_single_C = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C" and bond.order == Bond.SINGLE]
        # case A
        if len(attached_double_C) == 1 and len(attached_single_C) == 0:
            attached_double_C_O = [bond.a2 for bond in attached_double_C[0].bonds if bond.a2.elem == "O" and bond.a2.name != atom.name]
            if len(attached_double_C_O) == 1:
                attached_double_C_O_any = [bond.a2 for bond in attached_double_C_O[0].bonds if bond.a2.name != attached_double_C[0].name]
                if len(attached_double_C_O_any) == 0: return True
                else: return False
            else: return False
        # case B
        elif len(attached_double_C) == 0 and len(attached_single_C) == 1:
            attached_single_C_O = [bond.a2 for bond in attached_single_C[0].bonds if bond.a2.elem == "O" and bond.a2.name != atom.name]
            if len(attached_single_C_O) == 1:
                attached_single_C_O_any = [bond.a2 for bond in attached_single_C_O[0].bonds if bond.a2.name != attached_single_C[0].name]
                if len(attached_single_C_O_any) == 0: return True
                else: return False
            else: return False
        else: return False
    def is_charmm_OH1(atom):
        ''' For hydroxyl oxygen: oxygen bonded to a hydrogen and somthing else '''
        attached_H = [bond.a2 for bond in atom.bonds if bond.a2.elem == "H"]
        attached_other = [bond.a2 for bond in atom.bonds if not bond.a2.elem == "H"]
        if len(attached_H) == 1 and len(attached_other) >= 1: return True
        else: return False
    def is_charmm_OS(atom):
        ''' For ester oxygen: oxygen with 2 heavy atom single bonds, think this may work for ether too'''
        attached_H = [bond.a2 for bond in atom.bonds if bond.a2.elem == "H"]
        attached_other = [bond.a2 for bond in atom.bonds if not bond.a2.elem == "H"]
        if len(attached_H) == 0 and len(attached_other) >= 2: return True
        else: return False
    # ----- SULFER -----
    # ----- SULFER -----
    # ----- SULFER -----
    def is_charmm_S(atom):
        ''' For sulfur: sulfur with more than one bond '''
        num_bonds = len(atom.bonds)
        if num_bonds > 1: return True
        else: return False
    def is_charmm_SS(atom):
        ''' For thiolate sulfur: bare sulfur connected to a carbon '''
        num_bonds = len(atom.bonds)
        if num_bonds == 1: return True
        else: return False
    # main loop
    for a in atoms:
        if a.is_virtual: a.mm_type = "VIRT"
        elif a.elem == "H":
            assert(len(a.bonds) == 1) # hydrogen should only be attached to one atom
            at = a.bonds[0].a2
            if   is_charmm_H(a, at):   a.mm_type = "H"
            elif is_charmm_HC(a, at):  a.mm_type = "HC"
            elif is_charmm_HS(a, at):  a.mm_type = "HS"
            elif is_charmm_HB(a, at, peptoid):  a.mm_type = "HB"
            elif is_charmm_HA(a, at):  a.mm_type = "HA"
            elif is_charmm_HR1(a, at): a.mm_type = "HR1"
            elif is_charmm_HR2(a, at): a.mm_type = "HR2"
            elif is_charmm_HR3(a, at): a.mm_type = "HR3"
            elif is_charmm_HP(a, at):  a.mm_type = "HP"
            elif is_charmm_HE1(a, at): a.mm_type = "HE1"
            elif is_charmm_HE2(a, at): a.mm_type = "HE2"
            elif is_charmm_HF1(a, at): a.mm_type = "HF1"
            elif is_charmm_HF2(a, at): a.mm_type = "HF2"
            else: a.mm_type = "X"
        elif a.elem == "C" :
            if   is_charmm_CT(a):   a.mm_type = "CT"
            elif is_charmm_CT1(a):  a.mm_type = "CT1"
            elif is_charmm_CT2(a):  a.mm_type = "CT2"
            elif is_charmm_CT3(a):  a.mm_type = "CT3"
            elif is_charmm_CPH1(a): a.mm_type = "CPH1"
            elif is_charmm_CPH2(a): a.mm_type = "CPH2"
            elif is_charmm_CPT(a):  a.mm_type = "CPT"
            elif is_charmm_CY(a):   a.mm_type = "CY"
            elif is_charmm_CA(a):   a.mm_type = "CA"
            elif is_charmm_CC(a):   a.mm_type = "CC"
            elif is_charmm_CD(a):   a.mm_type = "CD"
            elif is_charmm_CS(a):   a.mm_type = "CS"
            elif is_charmm_CE1(a):  a.mm_type = "CE1"
            elif is_charmm_CE2(a):  a.mm_type = "CE2"
            elif is_charmm_CN(a):   a.mm_type = "CN"
            elif is_charmm_C(a):    a.mm_type = "C"
            else: a.mm_type = "X"
        elif a.elem == "N" :
            if   is_charmm_NY(a):  a.mm_type = "NY"
            elif is_charmm_NR1(a): a.mm_type = "NR1"
            elif is_charmm_NR2(a): a.mm_type = "NR2"
            elif is_charmm_NR3(a): a.mm_type = "NR3"
            elif is_charmm_N(a):   a.mm_type = "N"
            elif is_charmm_NH1(a): a.mm_type = "NH1"
            elif is_charmm_NH2(a): a.mm_type = "NH2"
            elif is_charmm_NH3(a): a.mm_type = "NH3"
            elif is_charmm_NC2(a): a.mm_type = "NC2"
            elif is_charmm_NP(a):  a.mm_type = "NP"
            elif is_charmm_NC(a):  a.mm_type = "NC"
            else: a.mm_type = "X"
        elif a.elem == "O" :
            if   is_charmm_O(a):   a.mm_type = "O" 
            elif is_charmm_OB(a):  a.mm_type = "OB" 
            elif is_charmm_OC(a):  a.mm_type = "OC" 
            elif is_charmm_OH1(a): a.mm_type = "OH1" 
            elif is_charmm_OS(a):  a.mm_type = "OS" 
            else: a.mm_type = "X"
        elif a.elem == "S" :
            if   is_charmm_SS(a):  a.mm_type = "SS"
            elif is_charmm_S(a):   a.mm_type = "S"
            else: a.mm_type = "X"
        elif a.elem == "P" : a.mm_type = "P"
        elif a.elem == "F" : a.mm_type = "X" # not doing halogens yet
        elif a.elem == "CL": a.mm_type = "X" # not doing halogens yet
        elif a.elem == "BR": a.mm_type = "X" # not doing halogens yet
        elif a.elem == "I" : a.mm_type = "X" # not doing halogens yet
        elif a.elem == "NA": a.mm_type = "NA"
        elif a.elem == "K" : a.mm_type = "K"
        elif a.elem == "MG": a.mm_type = "MG"
        elif a.elem == "FE": a.mm_type = "FE"
        elif a.elem == "CA": a.mm_type = "CA"
        elif a.elem == "ZN": a.mm_type = "ZN"
        else: a.mm_type = " X  " # this at least seems to be a legal MM atom type

def assign_partial_charges(atoms, net_charge=0.0):
    '''Assigns Rosetta standard partial charges, then
    corrects them so they sum to the desired net charge.
    Correction is distributed equally among all atoms.

    If non-zero partial charges were already assigned, no change is made.
    '''
    null_charge = [a for a in atoms if a.partial_charge is None]
    abs_charge = sum(abs(a.partial_charge) for a in atoms if a.partial_charge is not None)
    if len(null_charge) == 0 and abs_charge > 0:
        net_charge = sum(a.partial_charge for a in atoms if a.partial_charge is not None)
        print "Partial charges already fully assigned, no changes made; net charge %.3f" % net_charge
        return
    elif 0 < len(null_charge) and len(null_charge) < len(atoms):
        raise ValueError("Only some partial charges were assigned -- must be all or none.")
    std_charges = { # from Rosetta++ aaproperties_pack.cc
        "CNH2" : 0.550,
        "COO " : 0.620,
        "CH1 " : -0.090,
        "CH2 " : -0.180,
        "CH3 " : -0.270,
        "aroC" : -0.115,
        "Ntrp" : -0.610,
        "Nhis" : -0.530,
        "NH2O" : -0.470,
        "Nlys" : -0.620,
        "Narg" : -0.750,
        "Npro" : -0.370,
        "OH  " : -0.660,
        "Oaro" : -0.660, # copied from OH
        "ONH2" : -0.550,
        "OOC " : -0.760,
        "S   " : -0.160,
        "Nbb " : -0.470,
        "CAbb" : 0.070,
        "CObb" : 0.510,
        "OCbb" : -0.510,
        "Phos" : 1.500,
        "Hpol" : 0.430,
        "Hapo" : 0.095,
        "Haro" : 0.115,
        "HNbb" : 0.310,
        "H2O " : 0.000,
        "F   " : -0.250,
        "Cl  " : -0.130,
        "Br  " : -0.100,
        "I   " : -0.090,
        "Zn2p" : 2.000,
        "Fe2p" : 2.000,
        "Fe3p" : 3.000,
        "Mg2p" : 2.000,
        "Ca2p" : 2.000,
        "Na1p" : 1.000,
        "K1p " : 1.000,
        "VIRT" : 0.000,
    }
    curr_net_charge = 0.0
    for a in atoms:
        a.partial_charge = std_charges[ a.ros_type ]
        curr_net_charge += a.partial_charge
    # We only want to operate on non-virtual atoms now:
    atoms = [a for a in atoms if not a.is_virtual]
    charge_correction = (net_charge - curr_net_charge) / len(atoms)
    print "Total naive charge %.3f, desired charge %.3f, offsetting all atoms by %.3f" % (curr_net_charge, net_charge, charge_correction)
    curr_net_charge = 0.0
    for a in atoms:
        a.partial_charge += charge_correction
        curr_net_charge += a.partial_charge
    assert( abs(net_charge - curr_net_charge) < 1e-3 )

def assign_rotatable_bonds(bonds):
    '''Rotatable bonds are single bonds outside of rings
    with heavy atoms attached to both ends (i.e. not methyls) or non-C with one H.'''
    def is_ok(a):
        bonds_to_heavy = len(a.heavy_bonds)
        bonds_to_hydro = len([b for b in a.bonds if b.a2.is_H])
        return bonds_to_heavy >= 2 or (
            bonds_to_heavy == 1 and a.heavy_bonds[0].order == Bond.SINGLE
            and bonds_to_hydro == 1 and a.elem != "C"
        )
    for b in bonds:
        b.can_rotate = (
            b.order == Bond.SINGLE
            and not b.is_ring
            and is_ok(b.a1)
            and is_ok(b.a2)
        )
        b.is_proton_chi = (
            b.can_rotate
            and (len(b.a1.heavy_bonds) == 1 or len(b.a2.heavy_bonds) == 1)
        )
        b.mirror.can_rotate     = b.can_rotate
        b.mirror.is_proton_chi  = b.is_proton_chi

def assign_rigid_ids(atoms):
    '''Groups atoms that are connected in rigid units, i.e. no rotatable bonds.'''
    # Iterate through atoms, assigning them to rigids via depth first search
    def assign_to_rigid(atom, rig_id):
        if atom.rigid_id == rig_id: return # already visited this one
        assert(atom.rigid_id == 0)
        atom.rigid_id = rig_id
        # Recurse over all non-rotatable bonds and assign those atoms to same rigid
        for bond in atom.bonds:
            if not bond.can_rotate:
                assign_to_rigid(bond.a2, rig_id)
    num_rig_id = 0
    for atom in atoms:
        # We'll skip all but the first atom in each rigid as all atoms
        # in the rigid will be assigned by assign_to_rigid()
        if atom.rigid_id == 0:
            num_rig_id += 1
            assign_to_rigid(atom, num_rig_id)

def fragment_ligand(molfile):
    '''Sets Atom.fragment_id, Atom.conn_bonds, and Bond.connection_id.
    Returns number of fragments created, i.e. the largest valid fragment id.'''
    remaining_bonds = Set(molfile.bonds) # copy
    # Delete all split bonds from remaining_bonds
    # and number the connections they leave behind
    num_conn_id = 0
    for line in molfile.footer:
        if not line.startswith("M SPLT"): continue
        fields = line.split()
        atom1 = molfile.atoms[int(fields[2]) - 1]
        atom2 = molfile.atoms[int(fields[3]) - 1]
        bond_to_remove = None
        for b in remaining_bonds:
            if((b.a1 == atom1 and b.a2 == atom2)
            or (b.a1 == atom2 and b.a2 == atom1)):
                bond_to_remove = b
                break
        if b is None:
            raise ValueError("Cannot find bond to split between %s and %s" % (atom1.name, atom2.name))
        #elif b.can_rotate:
            #raise ValueError("Shouldn't split ROTATABLE bond between %s and %s" % (atom1.name, atom2.name))
        else:
            if b.can_rotate:
                print "WARNING: spliting ROTATABLE bond between %s and %s" % (atom1.name, atom2.name)
            print "Split bond between %s and %s" % (atom1.name, atom2.name)
            num_conn_id += 1
            bond_to_remove.connection_id = num_conn_id
            bond_to_remove.mirror.connection_id = bond_to_remove.connection_id
            bond_to_remove.a1.conn_bonds.append(bond_to_remove)
            bond_to_remove.a2.conn_bonds.append(bond_to_remove.mirror)
            remaining_bonds.remove(bond_to_remove)
    # Even though no single fragment can have more than 9 connections (CONN1 - CONN9),
    # the set of all fragments together may have more than 9 connections.
    #
    # Iterate through atoms, assigning them to fragments via depth first search
    def assign_to_fragment(atom, frag_id):
        if atom.fragment_id == frag_id: return # already visited this one
        assert(atom.fragment_id == 0)
        atom.fragment_id = frag_id
        # Recurse over all unbroken bonds and assign those atoms to same fragment
        for bond in atom.bonds:
            if bond in remaining_bonds or bond.mirror in remaining_bonds:
                assign_to_fragment(bond.a2, frag_id)
    num_frag_id = 0
    for atom in molfile.atoms:
        # We'll skip all but the first atom in each fragment as all atoms
        # in the fragment will be assigned by assign_to_fragment()
        if atom.fragment_id == 0:
            num_frag_id += 1
            assign_to_fragment(atom, num_frag_id)
    #for atom in molfile.atoms:
    #    print atom.name, atom.fragment_id
    # Assert that all atoms have been assigned to a fragment
    assert(len([a for a in molfile.atoms if a.fragment_id == 0]) == 0)
    # More than 9 fragments will break our current residue naming scheme,
    # which uses two letters plus a number from 1 to 9.
    if num_frag_id > 9:
        raise ValueError("More than 9 ligand fragments!")
    # No atom may have more than 1 connection outside the residue.
    # No two fragments may have more than one connection.
    frag_frag_conns = Set()
    for atom in molfile.atoms:
        if len(atom.conn_bonds) > 1:
            raise ValueError("Cannot create more than one connection at same atom (%s)" % atom.name)
        # Only check one permutation b/c the other gets caught going backwards along the same bond.
        for conn in atom.conn_bonds:
            pair1 = (conn.a1.fragment_id, conn.a2.fragment_id)
            #pair2 = (conn.a2.fragment_id, conn.a1.fragment_id)
            if pair1 in frag_frag_conns:
                #raise ValueError("Cannot create multiple connections between fragments %i and %i" % pair1)
                print "WARNING: Multiple connections between fragments (%i and %i) **NOT CURRENTLY SUPPORTED** by Rosetta!" % pair1
            frag_frag_conns.add(pair1)
            #frag_frag_conns.add(pair2)
    # Fragments should be about the same size as other residues
    for frag_id in range(1,num_frag_id+1):
        frag_atoms = [a for a in molfile.atoms if a.fragment_id == frag_id]
        num_atoms = len(frag_atoms)
        num_heavy_atoms = len([a for a in frag_atoms if not a.is_H])
        num_rot_bonds = len([b for b in molfile.bonds if b.a1.fragment_id == frag_id and b.a2.fragment_id == frag_id and b.can_rotate])
        if num_atoms < 3:
            # Mini-Rosetta atomtree requires at least 3 points to establish a coordinate system.
            print "Fragment %i: %s" % (frag_id, [a.name for a in frag_atoms])
            raise ValueError("Fragment %i has %i atoms; merge with another fragment or add virtual atoms to make 3 total" % (frag_id, num_atoms))
        if not (7 <= num_atoms <= 24):
            print "WARNING: fragment %i has %i total atoms including H; protein residues have 7 - 24 (DNA: 33)" % (frag_id, num_atoms)
        if not (4 <= num_heavy_atoms <= 22):
            print "WARNING: fragment %i has %i non-H atoms; protein residues have 4 - 14 (DNA: 22)" % (frag_id, num_heavy_atoms)
        if num_rot_bonds > 4:
            print "WARNING: fragment %i has %i rotatable bonds; protein residues have 0 - 4" % (frag_id, num_rot_bonds)
    print "Average %.1f atoms (%.1f non-H atoms) per fragment" % (
        float(len(molfile.atoms)) / float(num_frag_id), float(len([a for a in molfile.atoms if not a.is_H])) / float(num_frag_id))
    # Average stats tabulated by IWD from Richardson's Top500 database
    print "(Proteins average 15.5 atoms (7.8 non-H atoms) per residue)"
    return num_frag_id

def build_fragment_trees(molfile):
    '''Assigns a root atom for each fragment and parents and children.'''
    # Assign root atoms based on instructions in the molfile
    for line in molfile.footer:
        # Standard MDL style is with a space, but KWK has used "MROOT" in the past. Babel uses 2 spaces.
        if   line.startswith("M  ROOT"): molfile.atoms[ int(line.split()[2]) - 1 ].is_root = True
        elif line.startswith("M ROOT"):  molfile.atoms[ int(line.split()[2]) - 1 ].is_root = True
        elif line.startswith("MROOT"):   molfile.atoms[ int(line.split()[1]) - 1 ].is_root = True
    for frag_id in Set([a.fragment_id for a in molfile.atoms]):
        # If we want to have a default way of choosing the root atom, this is the place:
        root_atoms = [a for a in molfile.atoms if a.fragment_id == frag_id and a.is_root]
        if len(root_atoms) == 0:
            print "WARNING:  no root atom specified, using auto-selected NBR atom instead."
            (nbr, nbr_dist) = choose_neighbor_atom(molfile, frag_id)
            nbr.is_root = True
            root_atoms = [nbr]
        elif len(root_atoms) != 1:
            raise ValueError("You must have no more than one 'M ROOT' record in the molfile per ligand fragment")
        root_atom = root_atoms[0]
        # Protein residues appear to go depth first, so that all chi angles ride on each other.
        # Depth first assignment -- leads to very deep trees
        def tree_dfs(parent):
            # Want to visit non-H children first
            tmp_children = [b.a2 for b in parent.bonds]
            tmp_children.sort(lambda a,b: cmp(a.is_H, b.is_H))
            for child in tmp_children:
                if child.fragment_id != parent.fragment_id: continue
                if child.parent is not None or child.is_root: continue
                child.parent = parent
                parent.children.append(child)
                tree_dfs(child)
        tree_dfs(root_atom)
        ## Breadth first assigment -- minimizes depth of the tree
        #bfs_list = [root_atom]
        #while len(bfs_list) > 0:
        #    # pop first item
        #    parent = bfs_list[0]
        #    del bfs_list[0]
        #    for bond in parent.bonds:
        #        child = bond.a2
        #        if child.fragment_id != parent.fragment_id: continue
        #        if child.parent is not None or child.is_root: continue
        #        child.parent = parent
        #        parent.children.append(child)
        #        bfs_list.append(child)
        #    # sort heavy atom children before hydrogens
        #    parent.children.sort(lambda a,b: cmp(a.is_H, b.is_H))
    # Every atom should have a parent OR be a root
    assert(len([a for a in molfile.atoms if not a.is_root and a.parent is None]) == 0)

def assign_internal_coords(molfile):
    '''Sets up stubs/input_stubs and d,theta,phi for all atoms.'''
    for frag_id in Set([a.fragment_id for a in molfile.atoms]):
        root_atoms = [a for a in molfile.atoms if a.fragment_id == frag_id and a.is_root]
        assert(len(root_atoms) == 1)
        root_atom = root_atoms[0]
        # Assign stub atoms so we can calculate our internal coords
        # Root atom doesn't use input stubs -- has no internal coords
        def assign_stubs(me):
            # Find stub atoms for me.  Must store them b/c defined recursively.
            # Logic copied from core::kinematics::tree::Atom_, BondedAtom_, JumpAtom_ (*.hh)
            # Only the root atom is a jump in these cases, simplifying things somewhat.
            if me.is_root:
                me.stub1 = me
                me.stub2 = me.children[0] # first child
                if len(me.stub2.children) > 0:  me.stub3 = me.stub2.children[0] # first child of first child
                else:                           me.stub3 = me.children[1] # second child
            else:
                me.stub1 = me
                me.stub2 = me.parent
                me.stub3 = me.parent.stub2
                # Special case for first child of the root
                if me.parent.is_root and me.stub3 == me:
                    me.stub3 = me.parent.stub3
            #print "stubs", [x.name for x in (me, me.stub1, me.stub2, me.stub3)]
            # Assign input stubs to children and calculate internal coords
            prev_sibling = me.stub3
            parent = me # rename to make logic clearer
            # Have to store input_stub atoms b/c they're written to params file
            for child in parent.children:
                child.input_stub1 = parent.stub1
                child.input_stub2 = parent.stub2
                child.input_stub3 = prev_sibling # for first child, this is parent.stub3
                # Special case for second child of the root
                if parent.is_root and prev_sibling == parent.stub2:
                    #print "activate second child case! stub3 =", parent.stub3.name
                    child.input_stub3 = parent.stub3
                #print "input_stubs", [x.name for x in (child, child.parent, child.input_stub1, child.input_stub2, child.input_stub3)]
                # Now actually calculate spherical internal coordinates
                child.d, child.theta, child.phi = calc_internal_coords(child, child.input_stub1, child.input_stub2, child.input_stub3)
                # Recursive update of child's children
                assign_stubs(child)
                # Child is now previous sibling for next child in this loop
                prev_sibling = child
        # end assign_stubs()
        assign_stubs(root_atom)
        # Root has to have dummy input stub atoms to fill space in params file
        root_atom.input_stub1 = root_atom.stub1
        root_atom.input_stub2 = root_atom.stub2
        root_atom.input_stub3 = root_atom.stub3
        # Root has dummy values for d, theta, phi
        root_atom.d     = 0.0
        root_atom.theta = 0.0
        root_atom.phi   = 0.0
    # end loop over all fragments

def calc_internal_coords(child, input_stub1, input_stub2, input_stub3):
    '''Returns (d, theta, phi) for a point given it's three input stub points.'''
    #print "calc_internal_coords", [x.name for x in (child, input_stub1, input_stub2, input_stub3)]
    # Now actually calculate spherical internal coordinates
    # There's some weird "flip_stub" logic for first child of the root
    # iff theta is to be kept fixed (BondedAtom.cc) but I don't think it matters here.
    #
    # parent == parent.stub1 == child.input_stub1, always
    # (except for CONNECT atoms, where we use different stubs!!)
    d = r3.distance(child, input_stub1)
    if d < 1e-2: # very small d;  theta, phi don't make sense
        print "WARNING: very small d=%f for %s" % (d, child.name)
        theta = 0.0
        phi = 0.0
    else:
        theta = r3.angle( r3.from_to(input_stub2,input_stub1), r3.from_to(input_stub1,child) )
        if theta < 1e-2 or theta > 180 - 1e-2:
            # This always happens for first child of root:
            #print "WARNING: nearly parallel theta=%f for %s" % (theta, child.name)
            phi = 0.0
        else:
            phi = r3.dihedral(child, input_stub1, input_stub2, input_stub3)
    return (d, theta, phi)

def choose_neighbor_atom(molfile, frag_id):
    atoms = [atom for atom in molfile.atoms if atom.fragment_id == frag_id]
    bonds = [bond for bond in molfile.bonds if bond.a1.fragment_id == frag_id and bond.a2.fragment_id == frag_id]
    # If the entire fragment is rigid, max distance between any two atoms
    # is their Cartesian distance.
    # If not, an upper bound on the max distance between two atoms is the sum
    # of the bond lengths along the path between them that minimizes said sum.
    # We can make this bound tighter by doing Cartesian distances within a
    # rigid substructure and only summing bonds between rigid units.
    # We implement this by making all atoms within a rigid unit
    # pseudo-neighbors of each other, joined together by pseudo-bonds.
    #
    # Create an N x N array of all-against-all distances.
    # Syntax is awkward and Python lacks a reliable +Inf, so we use 1e100.
    na = len(atoms) # Number of Atoms
    all_all_dist = [ [1e100] * na for i in range(na) ]
    num_rot_bonds = len([b for b in bonds if b.can_rotate])
    if num_rot_bonds == 0: # easy case -- rigid fragment
        for i in range(0,na):
            all_all_dist[i][i] = 0
            for j in range(i+1,na):
                d = r3.distance( atoms[i], atoms[j] )
                all_all_dist[i][j] = d
                all_all_dist[j][i] = d
    else: # hard case -- flexible fragment
        # This will be queried several times so better to precompute it:
        nbrs = dict([ (a,Set()) for a in atoms ])
        for a in atoms:
            # a's neighbors:  bonded and in same fragment
            nbrs[a].update([b.a2 for b in a.bonds if b.a2.fragment_id == frag_id])
            # a's pseudo-neighbors:  in same rigid unit and same fragment
            nbrs[a].update([a2 for a2 in atoms if a2.rigid_id == a.rigid_id])
        for i in range(0,na):
            all_all_dist[i] = dijkstra(
                start = atoms[i],
                nodes = atoms,
                nbr = lambda a: nbrs[a],
                dist = r3.distance
            )
    # The "best" neighbor atom is the one with
    # the smallest {max distance to all other atoms}.
    max_dist = [ max(row) for row in all_all_dist ]
    best_idx = argmin(max_dist)
    return (atoms[best_idx], max_dist[best_idx])

def dijkstra(start, nodes, nbr, dist):
    '''Computes the shortest path from start node to all others.
    start - the node to start from
    nodes - the set of all nodes (includes start)
    nbr() - returns a list of a node's neighbors, nbr(n) = [n's neighbors]
    dist() - returns distance between two nodes that are neighbors, dist(n,m)
    Returns array of shortest distances to start in same order as input nodes.
    '''
    # Not as efficient as it could be, I think, but best I could manage.
    # Tmp objects [dist_from_start,node] sort properly
    DIST = 0; NODE = 1
    queue = [ [1e100,node] for node in nodes ] # 1e100  ~  +Inf
    # Allows lookup of best distance by name
    shortest = dict([ (q[NODE],q) for q in queue ])
    shortest[start][DIST] = 0 # start
    while len(queue) > 0:
        curr_idx = argmin(queue) # on first pass this is start
        curr = queue.pop(curr_idx)[NODE]
        curr_shortest = shortest[curr][DIST]
        for n in nbr(curr):
            new_dist = curr_shortest + dist(curr,n)
            if new_dist < shortest[n][DIST]:
                shortest[n][DIST] = new_dist
    return [ shortest[node][DIST] for node in nodes ]

def write_ligand_kinemage(f, molfile):
    if not isinstance(f, file): f = open(f, 'w')
    f.write("@text\n")
    f.write("View this file with KiNG or Mage from http://kinemage.biochem.duke.edu\n")
    f.write("@kinemage 1\n")
    f.write("@title {%s}\n" % os.path.basename(f.name))
    f.write("@onewidth\n")
    # Element markers
    elem_colors = { "C": "green", "N": "sky", "O": "red", "H": "gray", "S": "yellow", "P": "peach", "F": "bluetint", "CL": "cyan", "BR": "sea", "I": "lilac" }
    f.write("@balllist {element balls} color= gray radius= 0.1\n")
    for a in molfile.atoms: f.write("{%s} %s %.3f %.3f %.3f\n" % (a.elem, elem_colors.get(a.elem.upper(), "hotpink"), a.x, a.y, a.z))
    # All bonds, including which can rotate
    f.write("@vectorlist {all bonds} color= gray width= 1\n")
    for b in molfile.bonds:
        f.write("{%s}P %.3f %.3f %.3f\n" % (b.a1.name, b.a1.x, b.a1.y, b.a1.z))
        f.write("{%s}L %.3f %.3f %.3f\n" % (b.a2.name, b.a2.x, b.a2.y, b.a2.z))
    f.write("@vectorlist {rotatable bonds} color= white width= 4\n")
    for b in molfile.bonds:
        if not b.can_rotate: continue
        f.write("{%s}P %.3f %.3f %.3f\n" % (b.a1.name, b.a1.x, b.a1.y, b.a1.z))
        f.write("{%s}L %.3f %.3f %.3f\n" % (b.a2.name, b.a2.x, b.a2.y, b.a2.z))
    # Atom labels
    ai = index_atoms(molfile.atoms) # Atom order in the molfile
    f.write("@labellist {atom indices} color= white off\n")
    for a in molfile.atoms: f.write("{%i} %.3f %.3f %.3f\n" % (ai[a], a.x, a.y, a.z))
    f.write("@labellist {original names} color= white off\n")
    for a in molfile.atoms: f.write("{%s} %.3f %.3f %.3f\n" % (a.orig_name, a.x, a.y, a.z))
    f.write("@labellist {atom names} color= white off\n")
    for a in molfile.atoms: f.write("{%s} %.3f %.3f %.3f\n" % (a.name, a.x, a.y, a.z))
    f.write("@labellist {Rosetta types} color= white\n")
    for a in molfile.atoms: f.write("{%s} %.3f %.3f %.3f\n" % (a.ros_type, a.x, a.y, a.z))
    f.write("@labellist {MM types} color= white off\n")
    for a in molfile.atoms: f.write("{%s} %.3f %.3f %.3f\n" % (a.mm_type, a.x, a.y, a.z))
    f.write("@labellist {PDB names} color= white off\n")
    for a in molfile.atoms: f.write("{%s} %.3f %.3f %.3f\n" % (a.pdb_name, a.x, a.y, a.z))
    f.write("@labellist {partial charges} color= white off\n")
    for a in molfile.atoms: f.write("{%.2f} %.3f %.3f %.3f\n" % (a.partial_charge, a.x, a.y, a.z))
    # Each fragment and its atom tree, distinguished by color
    colors = ['deadwhite', 'purple', 'blue', 'sky', 'cyan', 'sea', 'green', 'lime', 'yellow', 'gold', 'orange', 'red']
    frag_ids = list(Set([a.fragment_id for a in molfile.atoms]))
    frag_ids.sort()
    for frag_id in frag_ids:
        color = colors[ frag_id % len(colors) ]
        f.write("@group {frag %i} master= {fragments}\n" % frag_id)
        # Find and mark root atom for this fragment (should only be one)
        f.write("@balllist {root atom} color= %s radius= 0.2 master={root atoms}\n" % color)
        root_atoms = [a for a in molfile.atoms if a.fragment_id == frag_id and a.is_root]
        for a in root_atoms: f.write("{%s frag %i} %.3f %.3f %.3f\n" % (a.name, a.fragment_id, a.x, a.y, a.z))
        (nbr, nbr_dist) = choose_neighbor_atom(molfile, frag_id)
        f.write("@ringlist {nbr atom} color= %s radius= %.3f master={nbr atoms} off\n" % (color, nbr_dist))
        f.write("{%s frag %i} %.3f %.3f %.3f\n" % (nbr.name, nbr.fragment_id, nbr.x, nbr.y, nbr.z))
        f.write("{%s frag %i} r=0.3 %.3f %.3f %.3f\n" % (nbr.name, nbr.fragment_id, nbr.x, nbr.y, nbr.z))
        f.write("@arrowlist {atom tree} color= %s master= {atom trees}\n" % color) # can specify width, angle, radius
        frag_atoms = [a for a in molfile.atoms if a.fragment_id == frag_id]
        for a in frag_atoms:
            if a.parent is None: continue
            a1 = a.parent; a2 = a
            f.write("{%s}P %.3f %.3f %.3f\n" % (a1.name, a1.x, a1.y, a1.z))
            f.write("{%s}L %.3f %.3f %.3f\n" % (a2.name, a2.x, a2.y, a2.z))
        # Connection arrows are half length and meet in middle of bond
        f.write("@arrowlist {connections} color= %s master= {connections}\n" % color) # can specify width, angle, radius
        for a1 in frag_atoms:
            for b in a1.conn_bonds:
                #a2 = r3.midpoint(b.a1, b.a2)
                # so they don't *quite* meet in the middle (45% across):
                a2 = r3.add( r3.mult(b.a1,0.55), r3.mult(b.a2,0.45) )
                f.write("{%s}P %.3f %.3f %.3f\n" % (a1.name, a1.x, a1.y, a1.z))
                f.write("{%s}L %.3f %.3f %.3f\n" % ("midpt", a2.x, a2.y, a2.z))
    f.close()

def write_param_file(f, molfile, name, frag_id, base_confs, max_confs):
    '''Writes a Molfile object to a file.
    f may be a file name or file handle.
    base_confs is the number of heavy-atom conformations generated by e.g. Omega
    max_confs is the maximum number of conformations desired after adding proton rotation
        The total number of confs may still be larger than max_confs,
        but at least we'll skip -ex# extra sampling of proton chis.
    '''
    close_file = False
    if not isinstance(f, file):
        f = open(f, 'w')
        close_file = True
    if frag_id == 1 and len(name) > 2: name = "%3.3s" % name
    else: name = "%2.2s%1i" % (name, frag_id)
    f.write("NAME %s\n" % name)
    f.write("IO_STRING %3.3s %1i\n" % (name, frag_id))
    f.write("TYPE LIGAND\n")
    f.write("AA UNK\n")
    atoms = [atom for atom in molfile.atoms if atom.fragment_id == frag_id]
    bonds = [bond for bond in molfile.bonds if bond.a1.fragment_id == frag_id and bond.a2.fragment_id == frag_id]
    root_atoms = [atom for atom in atoms if atom.is_root]
    assert(len(root_atoms) == 1)
    root_atom = root_atoms[0]
    # ORDER that atoms appear seems to imply the order of the atom tree,
    # rather than the order of the ICOOR_INTERNAL records (?)
    #for atom in atoms:
    #    f.write("ATOM %-4s %-4s %-4s %.2f\n" % (atom.name, atom.ros_type, atom.mm_type, atom.partial_charge))
    # So instead we write atoms depth-first, starting from root.
    def write_atoms(atom):
        f.write("ATOM %-4s %-4s %-4s %.2f\n" % (atom.name, atom.ros_type, atom.mm_type, atom.partial_charge))
        for a2 in atom.children: write_atoms(a2)
    write_atoms(root_atom)
    for bond in bonds:
        f.write("BOND %-4s %-4s\n" % (bond.a1.name, bond.a2.name))
    # Define chi angles
    # Iterating over the bonds is non-trivial and we need multiple passes, so we define a generator:
    def rot_bond_iter(bonds):
        '''Yields the tuples (bond, a, b, c, d) for each rotatable bond and the four chi atoms.'''
        for bond in bonds:
            if not bond.can_rotate: continue
            # define atoms a-b-c-d with d closest to tree root
            if bond.a1 == bond.a2.parent:
                b = bond.a2
                c = bond.a1
            elif bond.a2 == bond.a1.parent:
                b = bond.a1
                c = bond.a2
            else: raise ValueError("Rotatable bond %s - %s not in atom tree!" % (bond.a1.name, bond.a2.name))
            # Child MUST have a child for this to have counted as rotatable.
            # Might as well use the first one, it's more likely to be a heavy atom.
            # ... unless child is at end of tree and has no other children except
            #                                   |    |                  |
            # the CONN atom, e.g. an ether:  -> D -> C -> B(oxygen) -|> A(next frag)
            #                                   |    |                  |
            # This seems likely to blow Rosetta's mind, so let's not allow it:
            if len(b.children) < 1:
                raise ValueError("Can't define chi angle: no children for %s.  Don't split ether bonds, OK?  Try --no-param for debugging." % b.name)
            a = b.children[0]
            if c.parent is not None:
                d = c.parent
            else:
                # If c is root atom, have to choose a different anchor
                # c MUST have 2+ children for this bond to have been called rotatable
                # (if it doesn't, we made a logical error somewhere...)
                # So, use first child of c that's not b:
                # ... again, unless C is root and the O of an ether bond (no other kids)
                d = [k for k in c.children if k != b][0]
            yield (bond, a, b, c, d)
    def is_sp2_proton(a, b, c, d):
        '''Crude guestimate of H-bonding geometry'''
        # i.e., assume sp2 if the H's grandparent has anything other than single bonds to it.
        return ((a.is_H and len([bnd for bnd in c.heavy_bonds if bnd.order != Bond.SINGLE]) > 0)
             or (d.is_H and len([bnd for bnd in b.heavy_bonds if bnd.order != Bond.SINGLE]) > 0))
    # Do proton chi's first so we can use -ex1, -ex2, etc
    sorted_bonds = list(bonds) # make a copy then sort in place
    sorted_bonds.sort(key=lambda b: b.is_proton_chi, reverse=True)
    num_H_confs = base_confs
    for bond, a, b, c, d in rot_bond_iter(sorted_bonds):
        if bond.is_proton_chi:
            if is_sp2_proton(a, b, c, d): num_H_confs *= 6
            else: num_H_confs *= 9
    if num_H_confs > max_confs:
        print "WARNING: skipping extra samples for proton chis; would give %i conformers" % num_H_confs
    num_chis = 0
    for bond, a, b, c, d in rot_bond_iter(sorted_bonds):
        num_chis += 1
        # "Direction" of chi definition appears not to matter,
        # but we follow the convention of amino acids (root to tip)
        f.write("CHI %i %-4s %-4s %-4s %-4s\n" % (num_chis, d.name, c.name, b.name, a.name))
        if bond.is_proton_chi:
            # Only expand proton chis with extra samples if doing so won't generate tens of thousands of confs.
            if num_H_confs <= max_confs: extra = "1 20"
            else: extra = "0"
            if is_sp2_proton(a, b, c, d):
                f.write("PROTON_CHI %i SAMPLES 2 0 180 EXTRA %s\n" % (num_chis, extra))
            else:
                f.write("PROTON_CHI %i SAMPLES 3 60 -60 180 EXTRA %s\n" % (num_chis, extra))
    # Assign numbers to the connection points in this fragment
    num_conns = 0
    conn_nums = {}
    for a in atoms:
        for b in a.conn_bonds:
            num_conns += 1
            conn_nums[b] = num_conns
            if b.can_rotate: rot_flag = "CAN_ROTATE"
            else: rot_flag = "NO_ROTATE"
            f.write("CONNECT %-4s %s #CONN%i\n" % (a.name, rot_flag, num_conns));
    (nbr, nbr_dist) = choose_neighbor_atom(molfile, frag_id)
    f.write("NBR_ATOM %s\n" % nbr.name)
    f.write("NBR_RADIUS %f\n" % nbr_dist)
    # Convention seems to be a depth-first traversal from the root.
    # I don't know whether this matters, but it's the easy way anyhow.
    def write_icoords(a):
        f.write("ICOOR_INTERNAL   %-4s %11.6f %11.6f %11.6f  %-4s  %-4s  %-4s\n"
            % (a.name, a.phi, a.theta, a.d, a.input_stub1.name, a.input_stub2.name, a.input_stub3.name));
        # Have to re-create some of the internal coord logic here for connection points
        # The bad thing about the current arrangement is that CONN points
        # may be come after hydrogens and use them as reference points...
        prev_sibling = a.stub3
        for child in a.children:
            write_icoords(child)
            prev_sibling = child
        for conn in a.conn_bonds:
            inp_stub1 = a.stub1
            inp_stub2 = a.stub2
            inp_stub3 = prev_sibling
            # Special case when connection is 2nd child of root
            if inp_stub3 == inp_stub2:
                assert(a.is_root)
                inp_stub3 = a.stub3
            d, theta, phi = calc_internal_coords(conn.a2, inp_stub1, inp_stub2, inp_stub3)
            f.write("ICOOR_INTERNAL  CONN%1i %11.6f %11.6f %11.6f  %-4s  %-4s  %-4s\n"
                % (conn_nums[conn], phi, theta, d, inp_stub1.name, inp_stub2.name, inp_stub3.name));
            prev_sibling = conn.a2 # if we allowed multiple connects from one parent
    write_icoords(root_atom)
    #
    # XXX-FIXME:  still might need TYPE, PROPERTIES, FIRST_SIDECHAIN_ATOM, ACT_COORD_ATOMS
    #
    if close_file: f.close()

def write_poly_param_file(f, molfile, name, frag_id, peptoid):
    '''Writes a Molfile object to a file.
    f may be a file name or file handle.
    '''
    close_file = False
    if not isinstance(f, file):
        f = open(f, 'w')
        close_file = True
    else: name = "%2.2s" % (name)
    f.write("NAME %s\n" % name)
    f.write("IO_STRING %3.3s X\n" % (name))
    f.write("TYPE POLYMER\n")
    f.write("AA UNK\n")

    atoms = [atom for atom in molfile.atoms if atom.poly_ignore == False]
    bonds = [bond for bond in molfile.bonds if bond.a1.poly_ignore == False and bond.a2.poly_ignore == False]
    root_atoms = [atom for atom in atoms if atom.is_root]
    assert(len(root_atoms) == 1)
    root_atom = root_atoms[0]

    # write atoms
    for atom in atoms:
        if atom.poly_lower != True and atom.poly_upper != True:
            f.write("ATOM %-4s %-4s %-4s %.2f\n" % (atom.pdb_name, atom.ros_type, atom.mm_type, atom.partial_charge))
    # write bonds 
    for bond in bonds:
        if bond.a1.poly_lower != True and bond.a1.poly_upper != True and bond.a2.poly_lower != True and bond.a2.poly_upper != True:
            if bond.a1.elem == 'H':
                f.write("BOND %-4s %-4s\n" % (bond.a2.pdb_name, bond.a1.pdb_name))
            else:
                f.write("BOND %-4s %-4s\n" % (bond.a1.pdb_name, bond.a2.pdb_name))

    # write upper and lower connect
    f.write("LOWER_CONNECT N\nUPPER_CONNECT C\n")

    # Define chi angles
    # Iterating over the bonds is non-trivial and we need multiple passes, so we define a generator:
    def rot_bond_iter(bonds):
        '''Yields the tuples (bond, a, b, c, d) for each rotatable bond and the four chi atoms.'''
        for bond in bonds:
            if not bond.can_rotate: continue
            # define atoms a-b-c-d with d closest to tree root
            if bond.a1 == bond.a2.parent:
                b = bond.a2
                c = bond.a1
            elif bond.a2 == bond.a1.parent:
                b = bond.a1
                c = bond.a2
            else: raise ValueError("Rotatable bond %s - %s not in atom tree!" % (bond.a1.name, bond.a2.name))
            # Child MUST have a child for this to have counted as rotatable.
            # Might as well use the first one, it's more likely to be a heavy atom.
            # ... unless child is at end of tree and has no other children except
            #                                   |    |                  |
            # the CONN atom, e.g. an ether:  -> D -> C -> B(oxygen) -|> A(next frag)
            #                                   |    |                  |
            # This seems likely to blow Rosetta's mind, so let's not allow it:
            if len(b.children) < 1:
                raise ValueError("Can't define chi angle: no children for %s.  Don't split ether bonds, OK?  Try --no-param for debugging." % b.name)
            a = b.children[0]
            if c.parent is not None:
                d = c.parent
            else:
                # If c is root atom, have to choose a different anchor
                # c MUST have 2+ children for this bond to have been called rotatable
                # (if it doesn't, we made a logical error somewhere...)
                # So, use first child of c that's not b:
                # ... again, unless C is root and the O of an ether bond (no other kids)
                d = [k for k in c.children if k != b][0]
            yield (bond, a, b, c, d)
    def is_sp2_proton(a, b, c, d):
        '''Crude guestimate of H-bonding geometry'''
        # i.e., assume sp2 if the H's grandparent has anything other than single bonds to it.
        return ((a.is_H and len([bnd for bnd in c.heavy_bonds if bnd.order != Bond.SINGLE]) > 0)
             or (d.is_H and len([bnd for bnd in b.heavy_bonds if bnd.order != Bond.SINGLE]) > 0))

    sorted_bonds = list(bonds) # make a copy then sort in place
    sorted_bonds.sort(key=lambda b: b.is_proton_chi, reverse=False)

    # hacky check to see if chis in correct order else revese the order
    if not peptoid:
        for bond, a, b, c, d in rot_bond_iter(sorted_bonds):
            if a.poly_upper == False and b.poly_upper == False and c.poly_upper == False and d.poly_upper == False and a.poly_lower == False and b.poly_lower == False and c.poly_lower == False and d.poly_lower == False:
                print d.poly_n_bb, d.pdb_name, c.pdb_name, b.pdb_name, a.pdb_name
                if not d.poly_n_bb:
                    sorted_bonds.reverse()
                    print "REVERSED!!!"
                    break
                else:
                    break

    num_chis = 0
    for bond, a, b, c, d in rot_bond_iter(sorted_bonds):
        if a.poly_upper == False and b.poly_upper == False and c.poly_upper == False and d.poly_upper == False and \
                a.poly_lower == False and b.poly_lower == False and c.poly_lower == False and d.poly_lower == False:
            num_chis += 1
            # "Direction" of chi definition appears not to matter,
            # but we follow the convention of amino acids (root to tip)
            f.write("CHI %i %-4s %-4s %-4s %-4s\n" % (num_chis, d.pdb_name, c.pdb_name, b.pdb_name, a.pdb_name))
            if bond.is_proton_chi:
                # Only expand proton chis with extra samples if doing so won't generate tens of thousands of confs.
                extra = "1 20"
                if is_sp2_proton(a, b, c, d):
                    f.write("PROTON_CHI %i SAMPLES 2 0 180 EXTRA %s\n" % (num_chis, extra))
                else:
                    f.write("PROTON_CHI %i SAMPLES 3 60 -60 180 EXTRA %s\n" % (num_chis, extra))

    # choose neighbor atom first CB else CA
    nbr_list = [i for i,a in enumerate(molfile.atoms) if a.pdb_greek_dist == 'B']
    if len(nbr_list) >= 1:
        nbr_atom_index = nbr_list[0]
        nbr_atom = atoms[nbr_atom_index]
    else:
        for i,a in enumerate(molfile.atoms):
            if a.poly_ca_bb == True:
                nbr_atom_index = i
                nbr_atom = a
                break
    f.write("NBR_ATOM %s\n" % nbr_atom.pdb_name)

    # calc theoretical max neighbor radius
    na = len(molfile.atoms) # Number of Atoms
    all_all_dist = [ [1e100] * na for i in range(na) ]
    nbrs = dict([ (a,Set()) for a in molfile.atoms ])
    for a in molfile.atoms:
        nbrs[a].update([b.a2 for b in a.bonds])
    for i in range(0,na):
        all_all_dist[i] = dijkstra( start = molfile.atoms[i], nodes = molfile.atoms, nbr = lambda a: nbrs[a], dist = r3.distance )
    nbr_dist = max(all_all_dist[nbr_atom_index])
    f.write("NBR_RADIUS %f\n" % nbr_dist)

    #determin first side chain atom order of atoms should be n ca c o upper lower [side chain heavys] [hydrogens]
    non_bb_heavy_atoms = [a for a in atoms if a.poly_backbone == False and a.poly_lower == False and a.poly_upper == False and a.poly_ignore == False and a.elem != 'H']
    if len(non_bb_heavy_atoms) > 0:
        f.write("FIRST_SIDECHAIN_ATOM %s\n" % non_bb_heavy_atoms[0].pdb_name)
    else:
        f.write("FIRST_SIDECHAIN_ATOM NONE\n")

    # properties
    for line in molfile.footer:
        if line.startswith("M  POLY_PROPERTIES"):
            properties_list = line.split()[2:]
    f.write("PROPERTIES")
    for p in properties_list: f.write(" %s" % p)
    f.write("\n")

    # hacky hardcoding of stubs and iccords for backbone atoms and upper lower
    nbb = cabb = cbb = obb = upper = lower = 0
    for a in atoms:
        if a.poly_n_bb == True: nbb = a
        if a.poly_ca_bb == True: cabb = a
        if a.poly_c_bb == True: cbb = a
        if a.poly_o_bb == True: obb = a
        if a.poly_upper == True: upper = a
        if a.poly_lower == True: lower = a

    nbb.input_stub1,  nbb.input_stub2,  nbb.input_stub3  = nbb, cabb, cbb
    cabb.input_stub1, cabb.input_stub2, cabb.input_stub3 = nbb, cabb, cbb
    cbb.input_stub1,  cbb.input_stub2,  cbb.input_stub3  = cabb, nbb, cbb
    obb.input_stub1,  obb.input_stub2,  obb.input_stub3  = cbb, cabb, nbb
    upper.input_stub1, upper.input_stub2, upper.input_stub3 = cbb, cabb, nbb
    lower.input_stub1, lower.input_stub2, lower.input_stub3 = nbb, cabb, cbb

    nbb.d,  nbb.theta,  nbb.phi  = calc_internal_coords(nbb,  nbb.input_stub1,  nbb.input_stub2,  nbb.input_stub3)
    cabb.d, cabb.theta, cabb.phi = calc_internal_coords(cabb, cabb.input_stub1, cabb.input_stub2, cabb.input_stub3)
    cbb.d,  cbb.theta,  cbb.phi  = calc_internal_coords(cbb,  cbb.input_stub1,  cbb.input_stub2,  cbb.input_stub3)
    obb.d,  obb.theta,  obb.phi  = calc_internal_coords(obb,  obb.input_stub1,  obb.input_stub2,  obb.input_stub3)
    upper.d, upper.theta, upper.phi  = calc_internal_coords(upper, upper.input_stub1, upper.input_stub2, upper.input_stub3)
    lower.d, lower.theta, lower.phi  = calc_internal_coords(lower, lower.input_stub1, lower.input_stub2, lower.input_stub3)

    #write_icoords
    #for a in atoms:
    #    f.write("ICOOR_INTERNAL   %-4s %11.6f %11.6f %11.6f  %-4s  %-4s  %-4s\n" % (a.pdb_name, a.phi, a.theta, a.d, a.input_stub1.pdb_name, a.input_stub2.pdb_name, a.input_stub3.pdb_name));

    def write_icoords(a):
        if not a.poly_n_bb and not a.poly_ca_bb and not a.poly_c_bb:
            f.write("ICOOR_INTERNAL   %-4s %11.6f %11.6f %11.6f  %-4s  %-4s  %-4s\n" % (a.pdb_name, a.phi, a.theta, a.d, a.input_stub1.pdb_name, a.input_stub2.pdb_name, a.input_stub3.pdb_name));
        for child in a.children:
            if not child.poly_ignore:
                write_icoords(child)

    f.write("ICOOR_INTERNAL   %-4s %11.6f %11.6f %11.6f  %-4s  %-4s  %-4s\n" % (nbb.pdb_name, nbb.phi, nbb.theta, nbb.d, nbb.input_stub1.pdb_name, nbb.input_stub2.pdb_name, nbb.input_stub3.pdb_name));
    f.write("ICOOR_INTERNAL   %-4s %11.6f %11.6f %11.6f  %-4s  %-4s  %-4s\n" % (cabb.pdb_name, cabb.phi, cabb.theta, cabb.d, cabb.input_stub1.pdb_name, cabb.input_stub2.pdb_name, cabb.input_stub3.pdb_name));
    f.write("ICOOR_INTERNAL   %-4s %11.6f %11.6f %11.6f  %-4s  %-4s  %-4s\n" % (cbb.pdb_name, cbb.phi, cbb.theta, cbb.d, cbb.input_stub1.pdb_name, cbb.input_stub2.pdb_name, cbb.input_stub3.pdb_name));
    write_icoords(nbb)

    #
    # XXX-FIXME:  still might need TYPE, PROPERTIES, FIRST_SIDECHAIN_ATOM, ACT_COORD_ATOMS
    #
    if close_file: f.close()

def write_ligand_pdb(f, molfile_tmpl, molfile_xyz, resname, ctr=None, chain_id='X'):
    '''Writes a PDB file with a series of residues representing one ligand conformation.
    The topology (atom names, fragment divisions, etc) are taken from molfile_tmpl,
    while the actual XYZ coordinates are taken from molfile_xyz.
    resname provides the first two characters of the residue name.
    f may be a file name or file handle.'''
    if not isinstance(f, file): f = open(f, 'w')
    # If ctr is set, make it an offset vector for recentering ligands
    if ctr is not None:
        curr_ctr = r3.centroid([a for a in molfile_xyz.atoms if not a.is_H])
        ctr = r3.sub(ctr, curr_ctr)
    else: ctr = r3.Triple(0,0,0)
    atom_num = 0
    frag_ids = list(Set([a.fragment_id for a in molfile_tmpl.atoms]))
    frag_ids.sort()
    for frag_id in frag_ids:
        ai = index_atoms(molfile_tmpl.atoms) # 1-based index
        atoms = [a for a in molfile_tmpl.atoms if a.fragment_id == frag_id]
        for atom_tmpl in atoms:
            atom_xyz = molfile_xyz.atoms[ ai[atom_tmpl]-1 ]
            xyz = r3.add(atom_xyz, ctr)
            atom_num += 1
            if len(frag_ids) == 1 and len(resname) > 2:
                f.write("HETATM%5i %-4.4s %3.3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s  \n"
                    % (atom_num, atom_tmpl.name, resname,          chain_id, frag_id, xyz.x, xyz.y, xyz.z, 1.0, 20.0, atom_tmpl.elem))
            else:
                f.write("HETATM%5i %-4.4s %2.2s%1i %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s  \n"
                    % (atom_num, atom_tmpl.name, resname, frag_id, chain_id, frag_id, xyz.x, xyz.y, xyz.z, 1.0, 20.0, atom_tmpl.elem))
    f.write("TER"+(" "*77)+"\n")
    f.close()

def polymer_assign_backbone_atom_types(m):
    # first get POLY flags from molfile
    for line in m.footer:
        if line.startswith("M  POLY_N_BB"):
            m.atoms[ int(line.split()[2]) - 1 ].poly_backbone = True
            m.atoms[ int(line.split()[2]) - 1 ].poly_n_bb = True
        if line.startswith("M  POLY_CA_BB"):
            m.atoms[ int(line.split()[2]) - 1 ].poly_backbone = True
            m.atoms[ int(line.split()[2]) - 1 ].poly_ca_bb = True
        if line.startswith("M  POLY_O_BB"):
            m.atoms[ int(line.split()[2]) - 1 ].poly_backbone = True
            m.atoms[ int(line.split()[2]) - 1 ].poly_o_bb = True
        if line.startswith("M  POLY_C_BB"):
            m.atoms[ int(line.split()[2]) - 1 ].poly_backbone = True
            m.atoms[ int(line.split()[2]) - 1 ].poly_c_bb = True
        if line.startswith("M  POLY_UPPER"):
            m.atoms[ int(line.split()[2]) - 1 ].poly_backbone = True
            m.atoms[ int(line.split()[2]) - 1 ].poly_upper = True
        if line.startswith("M  POLY_LOWER"):
            m.atoms[ int(line.split()[2]) - 1 ].poly_backbone = True
            m.atoms[ int(line.split()[2]) - 1 ].poly_lower = True

def polymer_assign_backbone_atom_names(atoms, bonds, peptoid):
    ''' modifies the rosetta atom type of the hydorgen attached to the alpha carbon and backbone nitrogen to be the correct
    rosetta atom types '''
    # heavy atom atom names
    for atom in atoms:
        if atom.poly_ca_bb:
            atom.ros_type = "CAbb"
            atom.pdb_name = " CA "
        elif atom.poly_n_bb:
            atom.ros_type = "Nbb"
            atom.pdb_name = " N  "
            if peptoid:
                atom.mm_type = "NXX"
        elif atom.poly_o_bb:
            atom.ros_type = "OCbb"
            atom.pdb_name = " O  "
        elif atom.poly_c_bb:
            atom.ros_type = "CObb"
            atom.pdb_name = " C  "
        elif atom.poly_upper:
            atom.ros_type = "X"
            atom.pdb_name = "UPPER"
        elif atom.poly_lower:
            atom.ros_type = "X"
            atom.pdb_name = "LOWER"
    # alpha carbon hydrogen(s)
    for bond in bonds:
        if bond.a1.poly_ca_bb and bond.a2.is_H :
            bond.a2.ros_type = "Hapo"
            if not peptoid:
                bond.a2.pdb_name = " HA "
        elif bond.a1.is_H and bond.a2.poly_ca_bb :
            bond.a1.ros_type = "Hapo"
            if not peptoid:
                bond.a1.pdb_name = " HA "
    # backbone nitrogen hydrogen(s)
    for bond in bonds:
        if bond.a1.poly_n_bb and bond.a2.is_H :
            bond.a2.ros_type = "HNbb"
            bond.a2.pdb_name = " H  "
        elif bond.a1.is_H and bond.a2.poly_n_bb :
            bond.a1.ros_type = "HNbb"
            bond.a1.pdb_name = " H  "

def polymer_assign_ignored_atoms_bonds(m):
    ''' sets the ignore boolean for each atom in the list and for each bond with at least one atom in the list'''
    for line in m.footer:
        if line.startswith("M  POLY_IGNORE"):
            ignore_list = line.split()[2:]
    ignore_list = [int(i)-1 for i in ignore_list]
    #atoms
    for i,a in enumerate(m.atoms):
        if i in ignore_list:
            a.poly_ignore = True
    # bonds
    for bond in m.bonds:
        for i in ignore_list:
            if m.atoms[i] == bond.a1 or m.atoms[i] == bond.a2:
                bond.poly_ignore = True

def polymer_assign_pdb_like_atom_names_to_sidechain(atoms, bonds, peptoid):
    ''' Assign PDB like names to atoms. PDB names are based on the path distance from the alpha carbon
    Greek Alphabet: Alpha, Beta, Gamma, Delta, Epsilon, Zeta, Eta, Theta, Iota, Kappa, Lambda, Mu, Nu, Xi, Omicron, Pi, Rho, Sigma, Tau, Upsilon, Phi, Chi, Psi, Omega'''
    # greek alphabet eta, tao and omega are skipped because they are the same as previous letters
    greek_alphabet = ['A', 'B', 'G', 'D', 'E', 'Z', 'T', 'I', 'K', 'L', 'M', 'N', 'X', 'O', 'P', 'R', 'S', 'U', 'P', 'C']
    elem_atom_num = {'C': 6, 'N': 7, 'O': 8, 'F': 9, 'NA': 11, 'MG': 12, 'P':15, 'S':16, 'CL':17, 'K':19, 'CA':20, 'FE':26, 'ZN':30, 'BR':35, 'I':53}
    # find alpha carbon or the ""alpha nitrogen" for peptoids (still called ca_index below)
    print "PEPTOID" , peptoid
    if peptoid:
        for ca_index, atom in enumerate(atoms):
            if atom.poly_n_bb:
                print atom
                break
    else:
        for ca_index, atom in enumerate(atoms):
            if atom.poly_ca_bb:
                break
    print "CA or NA index is %d" % ca_index
    # assign heavy atom and hydrogen pdb_elem
    for atom in atoms:
        atom.pdb_elem = atom.elem
    # assign heavy atom pdb_greek_dist
    def path_dist(a,b):
        return 1
    na = len(atoms) # Number of Atoms
    all_all_dist = [ [1e100] * na for i in range(na) ]
    nbrs = dict([ (a,Set()) for a in atoms ])
    for a in atoms:
        #nbrs[a].update([b.a2 for b in a.bonds if b.a2.is_H == False and b.a2.poly_ignore == False and b.a2.poly_backbone == False])
        nbrs[a].update([b.a2 for b in a.bonds if b.a2.is_H == False and b.a2.poly_ignore == False and b.a2.poly_n_bb == False and b.a2.poly_c_bb == False and b.a2.poly_o_bb == False and b.a2.poly_upper == False and b.a2.poly_lower == False])
    for i in range(0,na):
        all_all_dist[i] = dijkstra( start = atoms[i], nodes = atoms, nbr = lambda a: nbrs[a], dist = path_dist )
    print "ALL TO ALL DIST"
    #debug
    for i,a in enumerate(atoms):
        print a, all_all_dist[i]
    print "ALL TO ALL DIST CA INDEX"
    print all_all_dist[ca_index] #DEBUG
    for i, a in enumerate(atoms):
        if peptoid:
            if not a.is_H and not a.poly_ignore and not a.poly_n_bb and not a.poly_c_bb and not a.poly_o_bb and not a.poly_upper and not a.poly_lower:
                print "ATOM: ", a
                print "DISTANCE: %f" % (all_all_dist[ca_index][i]-1) 
                a.pdb_greek_dist = greek_alphabet[all_all_dist[ca_index][i]-1]
        else:
            if not a.is_H and not a.poly_ignore and not a.poly_backbone:
                a.pdb_greek_dist = greek_alphabet[all_all_dist[ca_index][i]]
    debug = [a.pdb_greek_dist for a in atoms if not a.is_H ] #DEBUG
    print debug #DEBUG
    # assign heavy atom pdb_postfix_num (stupidly inefficient)
    def compare_atom_num(x,y):
        if elem_atom_num[atoms[x].elem] > elem_atom_num[atoms[y].elem]:
            return -1
        elif elem_atom_num[atoms[x].elem] == elem_atom_num[atoms[y].elem]:
            return 0
        elif elem_atom_num[atoms[x].elem] < elem_atom_num[atoms[y].elem]:
            return 1
    for i, g in enumerate(greek_alphabet):
        temp = [j for j,a in enumerate(atoms) if a.pdb_greek_dist == greek_alphabet[i]]
        if len(temp) > 1:
            temp.sort(compare_atom_num)
            for k, t in enumerate(temp):
                index = k+1
                atoms[t].pdb_postfix_num = "%d" % index
    #debug
    for a in atoms:
        if a.poly_ca_bb:
            print "DEBUG CA_BB: ", a, ":", a.pdb_prefix_num, ":", a.pdb_elem, ":", a.pdb_greek_dist, ":", a.pdb_postfix_num
            for b in a.bonds:
                print b
    # assign hydrogen pdb_greek_dist and pdb_postfix_num
    for a in atoms:
        if a.is_H and not a.poly_c_bb:
            a.pdb_greek_dist = a.bonds[0].a2.pdb_greek_dist
            a.pdb_postfix_num = a.bonds[0].a2.pdb_postfix_num
    # assign hydrogen pdb_prefix_num
    if peptoid:
        for a in atoms:
            if not a.is_H and not a.poly_ignore and not a.poly_n_bb and not a.poly_c_bb and not a.poly_o_bb and not a.poly_upper and not a.poly_lower:
                attached_h = [atoms.index(b.a2) for b in a.bonds if b.a2.is_H == True]
                for i,ah in enumerate(attached_h):
                    blah = i + 1
                    atoms[ah].pdb_prefix_num = "%d" % blah
    else:
        for a in atoms:
            if not a.is_H and not a.poly_backbone and not a.poly_ignore:
                attached_h = [atoms.index(b.a2) for b in a.bonds if b.a2.is_H == True]
                for i,ah in enumerate(attached_h):
                    blah = i + 1
                    atoms[ah].pdb_prefix_num = "%d" % blah
    # assign full pdb name
    for a in atoms:
        a.pdb_name = a.pdb_prefix_num + a.pdb_elem + a.pdb_greek_dist + a.pdb_postfix_num
    #debug
    for a in atoms:
        if a.poly_ca_bb:
            for b in a.bonds:
                print "DEBUG: ", b.a2.pdb_name

def polymer_reorder_atoms(molfile):
    ''' Reorders the atoms acording to the pdb ordering so that the order of the internal cords is correct '''
    def poly_atom_cmp(atom1, atom2):
        ''' Sorts based on special polymer backbone type, greek letter distance, postfix num, prefix num '''
        greek_alphabet = { ' ':0, 'A':1, 'B':2, 'G':3, 'D':4, 'E':5, 'Z':6, 'T':7, 'I':8, 'K':9, 'L':10, 'M':11, 'N':12, 'X':13, 'O':14, 'P':15, 'R':16, 'S':17, 'U':18, 'P':19, 'C':20, 'W':21, 'W':22, 'W':23, 'W':24, 'W':25, 'W':26, 'W':27, 'W':28, 'W':29, 'W':30, 'W':31 }
        # ignore
        if atom1.poly_ignore == True and atom2.poly_ignore == True: return 0
        elif atom1.poly_ignore == True and atom2.poly_ignore == False: return 1
        elif atom1.poly_ignore == False and atom2.poly_ignore == True: return -1
        elif atom1.poly_ignore == False and atom2.poly_ignore == False:
                    # special poly types
                    if atom1.poly_backbone == True and atom2.poly_backbone == False: return -1
                    elif atom1.poly_backbone == False and atom2.poly_backbone == True: return 1
                    elif atom1.poly_backbone == True and atom2.poly_backbone == True:
                        if atom1.poly_n_bb == True: return -1
                        elif atom1.poly_ca_bb == True and atom2.poly_n_bb != True: return -1
                        elif atom1.poly_c_bb == True and atom2.poly_n_bb != True and atom2.poly_ca_bb != True: return -1
                        elif atom1.poly_o_bb == True and atom2.poly_n_bb != True and atom2.poly_ca_bb != True and atom2.poly_c_bb != True: return -1
                        else: return 1
                    elif atom1.poly_backbone == False and atom2.poly_backbone == False:
                        # upper
                        if atom1.poly_upper == True and atom2.poly_upper == True: return 0 # this should never happen
                        if atom1.poly_upper == True and atom2.poly_upper == False: return 1
                        elif atom1.poly_upper == False and atom2.poly_upper == True: return -1
                        elif atom1.poly_upper == False and atom2.poly_upper == False:
                            # lower
                            if atom1.poly_lower == True and atom2.poly_lower == True: return 0 # this should never happen
                            if atom1.poly_lower == True and atom2.poly_lower == False: return 1
                            elif atom1.poly_lower == False and atom2.poly_lower == True: return -1
                            elif atom1.poly_lower == False and atom2.poly_lower == False:
                                # hydrogen
                                if atom1.elem == 'H' and atom2.elem != 'H': return 1
                                # greek distance
                                if greek_alphabet[atom1.pdb_greek_dist] < greek_alphabet[atom2.pdb_greek_dist]: return -1
                                elif greek_alphabet[atom1.pdb_greek_dist] > greek_alphabet[atom2.pdb_greek_dist]: return 1
                                elif greek_alphabet[atom1.pdb_greek_dist] == greek_alphabet[atom2.pdb_greek_dist]:
                                    # postfix num
                                    if atom1.pdb_postfix_num < atom2.pdb_postfix_num: return -1
                                    elif atom1.pdb_postfix_num > atom2.pdb_postfix_num: return 1
                                    elif atom1.pdb_postfix_num == atom2.pdb_postfix_num:
                                        # prefix num
                                        if atom1.pdb_prefix_num < atom2.pdb_prefix_num: return -1
                                        elif atom1.pdb_prefix_num > atom2.pdb_prefix_num: return 1
                                        else: return 0 #
    molfile.atoms.sort(poly_atom_cmp)
        

def main(argv):
    """
Converts a small molecule in an MDL Molfile with "M SPLT" and "M ROOT"
records into a series of .params residue definition files for Rosetta.
Also writes out the ligand conformation as PDB HETATMs.
If an SD file is given as input instead, the first entry is used for
generating topology / parameter files, and they all are used for
generating PDB-style coordinates in separate, numbered files.
Multiple models may also be supplied in MOL2 format, which does not support
M ROOT and M SPLT records but does allow for partial charges.
File type is deduced from the extension.

To divide a ligand into fragments by "breaking" bonds (optional):
M SPLT atom_no1 atom_no2

To specify a root atom for a ligand fragment (optional):
M ROOT atom_no

Expects that the input ligand has already had aromaticity "perceived",
i.e. that it contains aromatic bonds rather than alternating single and double
bonds (Kekule structure).

Optionally writes a kinemage graphics visualization of the atom tree,
neighbor atom selection, fragments, etc -- very helpful for debugging
and for visualizing exactly what was done to the ligand.
    """
    parser = OptionParser(usage="usage: %prog [flags] { INPUT.mol | INPUT.sdf | INPUT.mol2 }")
    parser.set_description(main.__doc__)
    # parser.add_option("-short", ["--long"],
    #   action="store|store_true|store_false",
    #   default=True|False|...
    #   type="string|int|float",
    #   dest="opt_name",
    #   help="store value in PLACE",
    #   metavar="PLACE",
    # )
    parser.add_option("-n", "--name",
        default="LG",
        help="name ligand residues NM1,NM2,... instead of LG1,LG2,...",
        metavar="NM"
    )
    parser.add_option("-p", "--pdb",
        default=None, # same as --name, see below
        help="prefix for PDB file names",
        metavar="FILE"
    )
    parser.add_option("-c", "--centroid",
        default=None, # same as --name, see below
        help="translate output PDB coords to have given heavy-atom centroid",
        metavar="X,Y,Z"
    )
    parser.add_option("-m", "--max-confs",
        default=5000, # 400 (default Omega max) * 9 (one sp3 H with -ex1) = 3600
        type="int",
        help="don't expand proton chis if above this many total confs",
        metavar="MAX"
    )
    parser.add_option("-k", "--kinemage",
        default=None,
        help="write ligand topology to FILE",
        metavar="FILE"
    )
    parser.add_option("--clobber",
        default=False,
        action="store_true",
        help="overwrite existing files"
    )
    parser.add_option("--no-param",
        default=False,
        action="store_true",
        help="skip writing .params files (for debugging)"
    )
    parser.add_option("--no-pdb",
        default=False,
        action="store_true",
        help="skip writing .pdb files (for debugging)"
    )
    parser.add_option("--polymer",
        default=False,
        action="store_true",
        help="write a polymer style param file instead of a ligand param file"
    )
    parser.add_option("--peptoid",
        default=False,
        action="store_true",
        help="modifier for the polymer flag, adjusts PDB style naming to be correct for peptoids"
    )
    (options, args) = parser.parse_args(args=argv)
    if options.pdb is None: options.pdb = options.name

    if len(args) < 1:
        parser.print_help()
        print "Must specify input .mol file!"
        return 1
    elif len(args) == 1:
        infile = args[0]
    else:
        parser.print_help()
        print "Too many arguments!"
        return 1

    ctr = None
    if options.centroid:
        f = options.centroid.split(",")
        if len(f) != 3:
            f = options.centroid.split()
            if len(f) != 3:
                print "Must say -centroid 'X,Y,Z'"
                return 5
        ctr = r3.Triple( float(f[0]), float(f[1]), float(f[2]) )
        #print "Centering ligands at %s" % ctr

    # There's a very strong order dependence to these function calls:
    # many depend on atom/bond variables set by earlier calls.
    infile_lc = infile.lower()
    if infile_lc.endswith(".mol2"):
        molfiles = read_tripos_mol2(infile)
    elif infile_lc.endswith(".mol") or infile_lc.endswith(".mdl") or infile_lc.endswith(".sdf"):
        molfiles = read_mdl_sdf(infile)
    else:
        print "Unrecognized file type, must be .mol/.sdf or .mol2!"
        return 6
    m = molfiles[0]
        # If -centroid not given, default is to center like first entry
    if ctr is None:
        ctr = r3.centroid([a for a in m.atoms if not a.is_H])
    print "Centering ligands at %s" % ctr
    add_fields_to_atoms(m.atoms)
    add_fields_to_bonds(m.bonds)
    find_virtual_atoms(m.atoms)
    uniquify_atom_names(m.atoms)
    check_bond_count(m.atoms)
    check_aromaticity(m.bonds)
    if options.polymer:
        polymer_assign_backbone_atom_types(m)
    assign_rosetta_types(m.atoms)
    assign_mm_types(m.atoms, options.peptoid)
    net_charge = 0.0
    for line in m.footer:
        # KWK's convention
        if line.startswith("M CHG"): net_charge = float(line[5:])
        # Official MDL format (integer only)
        # Oops, this isn't right!  MDL records N atom1 charge1 atom2 charge2 ...
        elif line.startswith("M  CHG"):
            charge_fields = line.split()[3:]
            net_charge = sum(int(c) for i,c in enumerate(charge_fields) if i%2 == 1)
        elif line.startswith("M  POLY_CHG"): net_charge = float(line.split()[2])
    assign_partial_charges(m.atoms, net_charge)
    assign_rotatable_bonds(m.bonds)
    assign_rigid_ids(m.atoms)
    num_frags = fragment_ligand(m)
    build_fragment_trees(m)
    assign_internal_coords(m)
    if options.polymer:
        print "Preforming polymer modifications"
        polymer_assign_ignored_atoms_bonds(m)
        polymer_assign_pdb_like_atom_names_to_sidechain( m.atoms, m.bonds, options.peptoid )
        polymer_assign_backbone_atom_names( m.atoms, m.bonds, options.peptoid )
        polymer_reorder_atoms(m)
    #uniquify_atom_names(m.atoms)
    if not options.no_param:
        for i in range(num_frags):
            if num_frags == 1: param_file = "%s.params" % options.pdb
            else: param_file = "%s%i.params" % (options.pdb, i+1)
            if not options.clobber and os.path.exists(param_file):
                print "File %s already exists -- aborting!" % param_file
                print "Use --clobber to overwrite existing files."
                return 2
            else:
                if options.polymer:
                    write_poly_param_file(param_file, m, options.name, 1, options.peptoid)
                    print "Wrote polymer params file %s" % param_file
                else:
                    write_param_file(param_file, m, options.name, i+1, len(molfiles), options.max_confs)
                    print "Wrote params file %s" % param_file
    if options.kinemage is not None:
        if not options.clobber and os.path.exists(options.kinemage):
            print "File %s already exists -- aborting!" % options.kinemage
            print "Use --clobber to overwrite existing files."
            return 3
        else:
            write_ligand_kinemage(options.kinemage, m)
            print "Wrote kinemage file %s" % options.kinemage
    if not options.no_pdb:
        for i, molfile in enumerate(molfiles):
            pdb_file = "%s_%04i.pdb" % (options.pdb, i+1)
            if not options.clobber and os.path.exists(pdb_file):
                print "File %s already exists -- aborting!" % pdb_file
                print "Use --clobber to overwrite existing files."
                return 4
            else:
                # m is used for names, molfile is used for XYZ
                write_ligand_pdb(pdb_file, m, molfile, options.name, ctr)
                print "Wrote PDB file %s" % pdb_file

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

    # Vestigal code for validating automatic atom typing:
    #m = read_mdl_molfile(sys.argv[1])
    #add_fields_to_atoms(m.atoms)
    #assign_rosetta_types(m.atoms)
    #for i, a in enumerate(m.atoms):
    #    err_flag = ""
    #    if a.name.strip() != a.ros_type.strip():
    #        if a.name == "CH1" and a.ros_type.startswith("CH"): pass
    #        else: err_flag = "********" #raise ValueError("typing mismatch!")
    #    print "%3i %4s --> %4s %s" % (i+1, a.name, a.ros_type, err_flag)

