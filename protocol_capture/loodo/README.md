LooDo (Loop-Directed Domain Insertion)
======================================

KEYWORDS: PROTEIN DESIGN

Author: Kristin M. Blacklock (kristin.blacklock@gmail.com)
Corresponding PI: Sagar D. Khare (khare@chem.rutgers.edu)

---

Publication describing the method:

* Blacklock KM, Yang L, Mulligan VK, Khare SD (2018) A computational method for the design of nested proteins by loop-directed domain insertion. 
Proteins. 86(3):354-369. doi: 10.1002/prot.25445.

## Brief
This protocol searches the conformational space for an insert domain using libraries of linker fragments for linker-to-parent domain superimposition followed by insert-to-linker superimposition. The relative positioning of the two domains (treated as rigid bodies) is sampled by a grid-based, mutual placement compatibility search.

## Executable ##
The LooDo application is implemented as a C++ executable in Rosetta.
   
`/path/to/Rosetta/source/bin/loodo.default.linuxgccrelease`


## Inputs ##
Three input files are required for the LooDo application:

1. __Inserted Domain__: A PDB file with no breaks in the chain, no "TER" lines 
in the file, and two residues on either end of the domain that will be used for
superimposition only (will not be present in the final structure).
2. __Parent Domain__: A PDB file with one break in its chain between the first 
and second halves of the domain where the inserted domain will be 
placed. Two residues on either side of the cutsite (total 4 residues) should be
present that will be used for superimposition (will not be present in the final
structure).
3. __Gridlig file__: A file describing a cloud of voxels in a user-defined
approximate position with respect to the parent domain that will be used to 
direct the placement of the inserted domain.

The inserted domain and parent domain PDBs can be generated in PyMOL, and
examples can be found in the inputs/ folder.

The gridlig file can be generated in the following way:

1. Open the parent domain in PyMOL. Fetch a PDB that is large enough to encompass the approximate area of the inserted domain.
  For example, type "fetch 1a3w" followed by "remove solvent" in PyMOL.

2. In PyMOL's editing mode, move the 1A3W structure 
to the approximate location where the inserted domain should be placed.

3. Save this new placement of the 1A3W object as a PDB file. This will be used as the "ligand" PDB in the gen\_lig\_grids application.

4. Run the gen\_lig\_grids application with your input files:
 
 `/path/to/Rosetta/source/bin/gen_lig_grids.default.linuxgccrelease -s <Parent Domain PDBfile> <"Ligand" PDBfile>`
        where the "Ligand" PDBfile is the 1A3W structure that has been moved near the parent domain.
 - Output should be a __ParentDomain.pdb_0.gridlig__ file.
 - There is a way to visualize this file in PyMOL for debugging purposes. 

## Running the Protocol ##

### Example Command Line ###
`/path/to/Rosetta/source/bin/loodo.default.linuxgccrelease @loodo.flags @general.flags`

### LooDo Flags ###

`-loodo:cap` (String) Cap or "Insert domain" PDB file

`-loodo:bot` (String) Bottom or "Parent domain" PDB file

`-loodo:fragAlength` (Integer) Length of the loop to insert between the N-terminal insertion site on the parent domain to the N-terminal insertion site on the insert domain, not including the two residues on either side of the insertion site that will be used for superimposition. Actual length will be `fragAlength + 4`. For example, setting this to `4` will mean an insertion of an 8-residue fragment, where the two flanking residues on either side of the insertion site will be replaced with two flanking residues of the loop fragment being inserted.

`-loodo:fragBlength` (Integer) Same as above but for the insertion site between the C-terminal insertion site on the insert domain and the C-terminal insertion site on the parent domain.

`-loodo:ins_begin` (Integer) This is the residue number (in pose numbering) of the last residue in the first half of the Parent Domain PDB, _including_ the two residues that will not be incorporated in the final structure. 

`-loodo:gridligpath` (String) /path/to/__ParentDomain.pdb_0.gridlig__

`-loodo:debug` (Bool) Setting this to `true` will dump all CapHit PDB files as they are produced by the protocol. It is advised to set this to `false`, cluster the CapHits by their Real6 values if there are too many to deal with in subsequent steps, and reconstitute the inserted domain placements using the `transform_loodo.default.linuxgccrelease` executable.

`-loodo:num_frags` (Integer) Sets the fragment picker to pick a certain number of loop fragments per frame, For example, if this flag is set to `4500` and `fragAlength` is set to `3`, the fragment picker would pick a total of `4500`*`7` fragments for fragment libraries A and B.

`-loodo:ca_ratio` (Real) Fraction of insert domain C-alphas required to occupy the active site grid.

`-loodo:com_in_grid` (Bool) If set to true, the C-alpha of the insert domain's center-of-mass residue is required to occupy the active site grid.

### General Flags ###
`-extra_res_fa` (String) path/to/params\_file for any ligands present in the input PDBs.

`-mute protocols.BumpGrid`

## Output ##

1. __CapHit_RT.txt__ file, which contains a list of descriptors for all passing insert domain placements.
2. __LinkerA*.pdb__ files, which are PDB files for each of the linkers from library A (N-terminal insertion site) that produced a passing insert domain placement.
3. __LinkerB*.pdb__ files, which are the same as above but for the C-terminal insertion site.
4. __Centered_Native.pdb__
5. If `-loodo:debug` is set to `true`, there may also be __CapHit*.pdb__ files.

## Next Steps

(Optional) If `-loodo:debug` was set to `false`, the Real6 values for each CapHit in the CapHit_RT.txt file can be clustered and reduced to a manageable number of hits.

Once a final set of CapHit RT's has been chosen,  use the `transform_loodo.default.linuxgccrelease` executable to reconstitute your clustered subset of CapHit RT's using the following command line:

`/path/to/Rosetta/source/bin/transform_loodo.default.linuxgccrelease -in:file:s NewCapHit_RT.txt -loodo:cap <Original Insert Domain PDBfile>`


#### Assemble Scaffolds for Downstream Protocols ####
Use the following python script in the `scripts/` directory to assemble the parent domain, LooDo-placed insert domains (CapHits), and linkers into a single PDB in correct residue order:

  `assemble_scaffolds.py -b <parent_domain_pdbfile> -i <insertion_site_resnum> -l <list_of_caphits_pdbfiles> -d <caphits_directory> -f <PyRosetta_init_flags>`
where the Parent Domain PDBfile is the same as the `-loodo:bot` input flag and the Insertion Site ResNum value is the same as the `-loodo:ins_begin` input flag. Only the CapHit PDB file names are required to be in the CapHits listfile. If the CapHits PDBs are in a different directory than where this script is run, use the `-d` flag to specify where the script should look for the CapHits. If there are any other flags you would like to supply to PyRosetta initialization (such as -ignore\_unrecognized\_res), use the `-f` flag to supply them (without the `-` character). 

The output will be of the form Scaffold\_CapHit\_*.pdb

Now the loops can be closed using your favorite loop closure protocol (the article linked above uses GenKic), followed by other optimization protocols.

** Note: Requires PyRosetta, which can be downloaded from http://www.pyrosetta.org