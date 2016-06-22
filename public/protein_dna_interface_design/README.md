/*!
The entire workflow for this demo should be described in a file
named README.dox.  It should describe an entire work flow, with
command lines, tested if possible.

@page protein_dna_interface_design Protein-DNA Interface Design (Demo)

@section documentation Documentation
See the rosettaDNA documentation at manual/applications/app_rosettaDNA.dox

@section pdb_input PDB Input Notes
Rosetta should be capable of reading in any protein-DNA complex in the modern PDB format
Some heavy-atom nucleotide residues may need to be converted by hand into canonical nucleotides. (e.g.: Iodocytosine...)
Alternative conformations will be ignored by Rosetta.

@section executable Rosetta Executable
bin/rosettaDNA.default[.{platform}{compiler}{mode}]

@section getting_started Getting Started
You need an XML-like RosettaScripts protocol file

@section scoring Scoring
[rosettaDNA.executable] -database [database path] -in:ignore_unrecognized_res -file:s starting_files/2h7h.pdb -score:weights dna -run:output_hbond_info -score:output_residue_energies -jd2:dd_parser -parser:protocol score.script -overwrite -out:prefix score_ > score.log

@section pack Packing

This packing routine finds optimal sidechain configurations for the wildtype amino acids in the vicinity of the target nucleotide base pair(s) (which are specified in the pack.script file as "dna_defs").

[rosettaDNA.executable] -database [database path] -in:ignore_unrecognized_res -file:s starting_files/2h7h.pdb -score:weights dna -run:output_hbond_info -score:output_residue_energies -jd2:dd_parser -parser:protocol pack.script -overwrite -out:prefix pack_ -ex1 > pack.log

@section design Design

This routine designs low energy amino acid identities and sidechains in the vicinity of the target nucleotides.

[rosettaDNA.executable] -database [database path] -in:ignore_unrecognized_res -file:s starting_files/2h7h.pdb -score:weights dna -run:output_hbond_info -score:output_residue_energies -jd2:dd_parser -parser:protocol design.script -overwrite -out:prefix design_ -ex1 > design.log

@section design Multi-state Design

This routine optimizes a starting population of different single-state designs for specificity toward the target DNA sequence vs. its competitors, using a genetic algorithm that involves mutation and recombination of the most sequence-specific solutions over multiple generations.

[rosettaDNA.executable] -database [database path] -in:ignore_unrecognized_res -file:s starting_files/2h7h.pdb -score:weights dna -run:output_hbond_info -score:output_residue_energies -jd2:dd_parser -parser:protocol multistate.script -overwrite -out:prefix multistate_ -ex1 > multistate.log

Typically, this multi-state genetic algorithm routine should be run with a population >= 1000 protein sequences, for >= 50 generations.

*/
