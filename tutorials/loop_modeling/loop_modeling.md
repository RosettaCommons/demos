Loop Modeling
=============
KEYWORDS: LOOPS HOMOLOGY_MODELING GENERAL   
Tutorial by Shourya S. Roy Burman (ssrb@jhu.edu)    
Created 23 June 2016

[[_TOC_]]

Summary
-------
Rosetta has a quite a few methods for modeling peptide segments in proteins. By the end of this tutorial, you should be able to understand:

* The variety of loop modeling methods in Rosetta
* How to join a break in the protein chain (without missing residues)
* How to model missing segments in proteins
* How to extend or delete peptide segments at the termini
* How to refine segments in proteins
* How to combine loop modeling with other protocols
* How to use loop modeling on non-canoncials

>This tutorial will not cover the algorithmic details of the loop modeling methods. You will be directed the documentation explaining the algorithm.

Loop Modeling Methods
---------------------
>Loop Modeling in not restricted to segments with a blank DSSP secondary structure assignment. They are more generally applicable to any short peptide fragment joining larger protein segments.

* CCD
* KIC
* NGK
* Generalized KIC
* Remodel
* Loop Hash

Closing Breaks in Protein Chains
--------------------------------
Sometimes you have chain breaks in your protein, perhaps when threading from a homolog. In this case, there are no missing residues, just that the backbone itself is not closed. An example input PDB is provided at `<path_to_Rosetta_directory>/demos/tutorials/loop_modeling/input_files/3gbn_Ab.pdb` where the connection between residue numbers 127 and 128 is severed. To fix this, we will use the kinematic closure protocol. First, we need to write a short _loop file_, detailing which residues as to be modeled as loops and where the cutpoint is.

    LOOP 125 130 0 0 1
    
This 2<sup>nd</sup> and the 3<sup>rd</sup> columns in this file correspond the the loop start and end residues. The 4<sup>th</sup> column indicates the cut point where the chain is broken. It must be between (including) the start and end residues. 0 is the default option which allows Rosetta to pick a cut point. Note that **the residue numbering in the loop file is not based on PDB numbering but on Rosetta internal numbering.**. The 4th column represents the skip rate, which we have set to 0. Setting the last column to 1 makes Rosetta start building from an extended structure.

To close the loop, run:

    $> <path_to_Rosetta_directory>/main/source/bin/loopmodel.linuxgccrelease @flag_basic_KIC 

    
This should take about 2 minutes at the end of which you will produce one PDB and one score file in `output_files`. This PDB will now have a closed loop.

A list of further options and documentation can be found [here](https://www.rosettacommons.org/docs/latest/application_documentation/structure_prediction/loop_modeling/loopmodel-kinematic).

Modeling Missing Loops
----------------------
Modeling missing loops is a difficult problem. Rosetta has a variety of algorithms for this purpose.


Extending and Deleting the Termini
----------------------------------

Refining Peptide Segments
-------------------------

Combining Loop Modeling with other Protocols
--------------------------------------------

```
<ROSETTASCRIPTS>
    <TASKOPERATIONS>
    </TASKOPERATIONS>
    <SCOREFXNS>
        <tala weights=talaris2014 symmetric=0 />
    </SCOREFXNS>
    <FILTERS>
    </FILTERS>
    <MOVERS>
        <PeptideStubMover name=pep_stub reset=1>
            <Append resname="CYD" />
            <Append resname="ALA" />
            <Append resname="ALA" />
            <Append resname="ALA" />
            <Append resname="ALA" />
            <Append resname="ALA" />
            <Append resname="ALA" />
            <Append resname="ALA" />
            <Append resname="ALA" />
            <Append resname="CYD" />
        </PeptideStubMover>
        
        <SetTorsion name=tor1>
            <Torsion residue=ALL torsion_name=phi angle=-64.7/>
            <Torsion residue=ALL torsion_name=psi angle=-41.0/>
            <Torsion residue=ALL torsion_name=omega angle=180/>
        </SetTorsion>

	<SetTorsion name=tor2>
            <Torsion residue=pick_atoms angle=random>
                <Atom1 residue=1 atom="N"/>
                <Atom2 residue=1 atom="CA"/>
                <Atom3 residue=1 atom="CB"/>
                <Atom4 residue=1 atom="SG"/>
            </Torsion>
            <Torsion residue=pick_atoms angle=random>
                <Atom1 residue=10 atom="N"/>
                <Atom2 residue=10 atom="CA"/>
                <Atom3 residue=10 atom="CB"/>
                <Atom4 residue=10 atom="SG"/>
            </Torsion>
	</SetTorsion>
        
        <DeclareBond name=bond res1=1 atom1="SG" res2=10 atom2="SG"/>

	<GeneralizedKIC name=genkic closure_attempts=200 stop_when_n_solutions_found=0 selector="random_selector">
		<AddResidue res_index=3 />
                <AddResidue res_index=2 />
                <AddResidue res_index=1 />
                <AddResidue res_index=10 />
                <AddResidue res_index=9 />
                <AddResidue res_index=8 />
                <AddResidue res_index=7 />
                <AddResidue res_index=6 />
                <AddResidue res_index=5 />
		<SetPivots res1=3 atom1="CA" res2=1 atom2="SG" res3=5 atom3="CA" />
		<CloseBond prioratom_res=10 prioratom="CB" res1=10 atom1="SG" res2=1 atom2="SG" followingatom_res=1 followingatom="CB" bondlength=2.05 angle1=103 angle2=103 randomize_flanking_torsions=true />
		<AddPerturber effect="randomize_alpha_backbone_by_rama">
			<AddResidue index=3 />
			<AddResidue index=2 />
			<AddResidue index=1 />
			<AddResidue index=10 />
			<AddResidue index=9 />
            <AddResidue index=8 />
            <AddResidue index=7 />
            <AddResidue index=6 />
            <AddResidue index=5 />
		</AddPerturber>
		<AddFilter type="loop_bump_check" />
	</GeneralizedKIC>
        
    </MOVERS>
    <APPLY_TO_POSE>
    </APPLY_TO_POSE>
    <PROTOCOLS>
        <Add mover=pep_stub/>
        <Add mover=tor1/>
        <Add mover=tor2/>
        <Add mover=bond/>
        <Add mover=genkic/>
    </PROTOCOLS>
</ROSETTASCRIPTS>
```

