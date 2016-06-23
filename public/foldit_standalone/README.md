Foldit Standalone Protocol Capture
==================================

KEYWORDS: UTILITIES GENERAL

Authors: Seth Cooper, David Kim, Firas Khatib, Adrien Treuille (PI), Janos 
Barbero, Jeehyung Lee, Michael Beenen, Andrew Leaver-Fay, Matthew Smith 
(Speaker), David Baker (PI), Zoran Popovic (PI)

Title of Talk: Standalone Foldit: A New and Powerful Interface to Rosetta

Date of Talk: Wednesday, August 4, 2010

Session: Interfaces to Rosetta

---

This protocol describes how to obtain the binaries and code for Foldit 
Standalone.  It also includes the example files used to demonstrate the 
software as well as the scripts and wrappers used to make the scripts.

Protein design with Rosetta typically involves transitioning between
command line tools (Rosetta) to generate designs, and graphical tools
(PyMol) to evaluate designs.  By tightening the feedback loop between
design and evaluation, we may be able to improve our design
methodology.  Standalone Foldit allows the protein designer to make
mutations; minimize ligands, sidechains, and the protein back bone;
and remodel loops, all with a sophisticated scripting interface,
advanced undo functionality, and the Rosetta energy function.

Additionally, Foldit can be used to evaluate and modify protein, DNA,
and RNA structures and any mixture of these.  Standalone Foldit is
under active development and has been used to design active enzymes.
This new interface to Rosetta is a great balance between usability,
interactivity, and power.

Running the protocol
--------------------

FolditStandalone is not a command line protocol, it is a graphical interface to 
Rosetta.

Versions
--------

Standalone Foldit demo was done with svn 35710, newer revisions should work 
similarly.

References
----------

Cooper S, Khatib F, Treuille A, Barbero J, Lee J, Beenen M, Leaver-Fay A, Baker 
D, Popovic Z, Foldit players. Predicting protein structures with a multiplayer 
online game. Nature 2010;466(7307):756-760.

Other comments
--------------

To use Lua script for reverting back to native sequence:

Place the attached scripts somewhere in your path (e.g. your scripts 
directory), and make all four executable by typing, in a terminal:

    chmod u+x lua_to_orig revert_to_orig.lua pos_id_lua.py pdb_fasta.pl

To generate the reversion script:

    lua_to_orig scaffold.pdb > /absolute/lua_script_location/lua_script_name

Then, at the lua prompt in Foldit, type:

    dofile("/absolute/lua_script_location/lua_script_name")
