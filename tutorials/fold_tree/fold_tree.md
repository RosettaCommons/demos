# Fold tree

In this tutorial you will learn about the concept of fold tree, internal coordinates, and how to use them wisely to obtain meaningfull outputs.

#### Internal Coordinates
If you open a PDB file and look at it, you can see that it has the information of all the atoms with their (x,y,z) coodrinates in the 3D space. When you want to move a protein, you can change all 3 coordinates of each atom. In other words, each atom has 3 degrees of freedom. However, 3 degrees of freedom is simply too much for our calculations. In order to make things simpler, Rosetta uses **internal coordinates** instead. You can imagin that another way of defining an atom uniquely is to use the bond, angle, and torsion of it with regards to its neighbor atoms. When you use bond, angle, and torsion angle values of residues, instead of the x,y,z values, to describe an atom, you are using its **internal coordinates**. 

In the internal coordinate world, you move objects by changing the bond, angle, and torsions, so there are again 3 degrees of freedom. However, within proteins and for most of the applications we are interested in, the bonds an the angles between atoms remain unchanged. So, the motions are mostly happening by a change in the torsion angles. Hence, we are reducing the degrees of freedom to 1.

If you go to
```
<path-to-Rosetta>/Rosetta/main/database/chemical/residue_type_sets/fa_standard/residue_types/l-caa
```
you can find several [params](https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/Residue-Params-file) files for all 20 amino acids. Open one and look at it. The lines started with ICOOR_INTERNAL show the internal coordinates of the atoms in that residue (see [here](https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/Residue-Params-file) for details of what each column in ICOOR_INERNAL line means.

#### The Lever Arm Effect
When you change the torsion angles in a given sets of residues but you set the rest to not move, what will happen? Let's see the simple scheme below:

![Figure for changes in the torsion](https://github.com/RosettaCommons/demos/blob/hssnzdh2/parisa_XRW/tutorials/figures/small_moves.png)

You can see that in order to change one specific torsion angle in the cyan residue, all of the atoms that are to the right of that change should all move. In other words, although they are held fixed "WITH RESPECT TO EACH OTHER", they still move together to accommodate changes in the _downstream_ cyan residue.

Now let's look at the same example, but with a small difference:

![figure for bigger protein changes](https://github.com/RosettaCommons/demos/blob/hssnzdh2/parisa_XRW/tutorials/figures/big_moves.png)

You can see that a small change in one torsion angle at the N-terminal of the protein can cause huge movements in the rest of the protein. This is called **the lever arm effect**.
How can we avoid this? 

#### Using the Fold Tree
In this part we see how to use *fold tree* to control movements of the parts in a protein with respects to each other. The **fold tree** is a way to tell Rosetta the connectivity between residues in a given structure. It defines what residues are _UPSTREAM_ or _children_ and what residues are _DOWNSTREAM_ or _parents_. Let's go through some examples to make things clear.

In this example, you will be using the structure of the N-terminal part of nucleocapsi from coronaviruse. Navigate to the fold-tree directory in tutorials using this command:
```
> cd <path-to-Rosetta>/Rosetta/demos/tutorials/fold_tree
```
Open the capsid.pdb file provided for you in the inputs directory. You can see that it has a long, unstructured N-terminal and the rest of the structure is mostly beta sheets. We know that the N-terminal region is disordered, so before working with the PDB, we want to relax just that part.
The cps_relax1.xml is a [Rosetta script] provided to you that relaxes residues 1-20 of the structure using a [MoveMap]. It also uses the default Rosetta fold tree, which is shown in caps_tree1.ft file in the input directory: 
```
FOLD_TREE EDGE 1 133 -1
```
You start from residue 1 in chain A and go all the way to the last residue (133 in this example) of the chain. So any movements in the downstream residues (residues that are close to N-terminal) will cause movements in residues upstream (closer to C-terminal). This chain is a non-disrupted continuous cluster of residues that are connected by covalent bonds in a linear way. We call this cluster an **EDGE**. **-1** means that the residues in the EDGE are connected through covalent bonds. 

In your fold_tree folder, run the script using the command below:
```
$> ../../../main/source/bin/rosetta_scripts.default.linuxgccrelease -in:file:s inputs/capsid.pdb -parser:protocol inputs/caps_relax1.xml -out:prefix test1_
```
NOTE: you may need to change your executable from what is provided here depending on your compilation. Check [here].

After 5-10 minutes, you can see that the output (test1_capsid_0001.pdb) is generated. 
Now compare the original structure with the output. You can see that the N-terminal part has moved drastically compared to the original structure. If we align the N-terminal parts, you can see that the rest of the protein has moved drastically to accomodate changes in the small 20 residue N-terminal.

![showing the monomer movement](https://github.com/RosettaCommons/demos/blob/hssnzdh2/parisa_XRW/tutorials/figures/ubq1_test1.png)

Now, let's see what happens when we sweep the position of the downstream and upstream, telling Rosetta that movements should propagate from C to N terminal. This is shown in caps_tree2.ft file:
```
FOLD_TREE EDGE 133 1 -1
```
You can see that the fold tree is directional.

let's run one round of relax on the capsid.pdb structure but this time with the new fold tree:
```
$> ../../../main/source/bin/rosetta_scripts.default.linuxgccrelease -in:file:s inputs/capsid.pdb -parser:protocol inputs/caps_relax1.xml -out:prefix test1_
```
This shouldn't take longer than 10 min. Look at the test2_capsid_0001.pdb output structure and compare it with both the capsid.pdb and the test1_capsid_0001.pdb structure. (the two outpus are provided in the outputs directory). You can see now that the movements in the protein are much slighter. This is because now the C-terminal parts are downstream of where the changes are happening and are not moving as a result of movements in the N-terminal.

![showing monomer improved](https://github.com/RosettaCommons/demos/blob/hssnzdh2/parisa_XRW/tutorials/figures/ubq_test2.png)

So if you have a protein with very flexible N-terminal that will move a lot during your run, you may want to change the fold tree.

#### Complexes With More Than One Chain and Jumps

Now let's run another example. In this example we are running a hypothetical ubiquitin dimer. The pdb is called ubq_dimer.pdb in the inputs directory.

Let's see what happens if I relax the pose. In the ubq_relax1.xml [Rosetta script] I used a [movemap] to fix the backbone expect for the small loop in the C-terminal. In the first run, we will use the default fold tree of Rosetta for protein complexes. Let's take a look at ubq_tree1.ft that contains this fold tree:
```
FOLD_TREE EDGE 1 76 -1 EDGE 76 77 1 EDGE 77 152 -1
```
You can see that we have three **EDGEs**. The first one should be familiar. It says the first edge start from residue 1 and goes through residu 76 by covalent linkage (basically, chain A). The third one is doing the same thing but for the chain B.(Note that these are pose numberings. See [here]). Now let's look at the second edge. You can see that it goes from the last residue of chain A (76) to the first residue of chain B (77). This is not a covalent interaction but we want Rosetta to know these are part of the same complex and they should remain somewhat close. So, we connect this two with an imaginary link that we call a **JUMP**. Each jump is an edge by tiself. Jumps are shown by positive numbers. The first jump is shown by 1 and the numbers go up for each jump.

Now let's relax our dimer using this default fold tree. Run this command in your fold-tree directory:

```
$> ../../../main/source/bin/rosetta_scripts.default.linuxgccrelease -s inputs/ubq_dimer.pdb -parser:protocol inputs/ubq_relax1.xml -out:prefix test1_
```
After 5-10 min, you should have an output called test1_ubq_dimer_0001.pdb. Take a look at the structure and compare it to the original dimer. You can see that the second chain moved a lot from it's original position. This is another example of the **lever effect**. The movements in the last few residues in C_terminal of chain A has propagated all the way to chain A.

![dimer with default fold_tree](https://github.com/RosettaCommons/demos/blob/hssnzdh2/parisa_XRW/tutorials/figures/ubq_dimer1.png)

Now, let's see how we can fix it. I know based on my experiments that the Val71 in chain A and Val146 in chain B are important in the interface formation and I want to make sure they stay close during relaxation. So, I want to re-define my fold tree so that residues 71 and 146 are the immediate parents, or the most downstream. Take a look at the modified fold tree (inputs/ubq_tree2.ft)
```
FOLD_TREE EDGE 71 1 -1 EDGE 71 76 -1 EDGE 71 146 1 EDGE 146 77 -1 EDGE 146 152 -1
```
You can see that now I have multiple **EDGEs**. Chain A is now defined by two edge: the first one goes through from residue 71 all the way back to N-terminal of chain A and the second 1 goes from residue 71 to C-terminal of chain A. The **JUMP** is now from 71 to 146, the two parents. And then the chain B is also defined by two edges now. The scheme below shows how the two fold tree differ:

![scheme of fold tree](https://github.com/RosettaCommons/demos/blob/hssnzdh2/parisa_XRW/tutorials/figures/FT_scheme.png)

Now, let's run again, but this time with the new fold tree we just defined:
```
$> ../../../main/source/bin/rosetta_scripts.default.linuxgccrelease -s inputs/ubq_dimer.pdb -parser:protocol inputs/ubq_relax2.xml -out:prefix test2_
```
Take a look at the test2_ubq_dimer_0001.pdb output and compare it with the original one. You can see that the drastic movements of chain B due C-terminal relaxation of chain A is now diminished.

![dimer with modified fold tree](https://github.com/RosettaCommons/demos/blob/hssnzdh2/parisa_XRW/tutorials/figures/ubq_dimer2.png)

#### Final Points
Now you know what a fold tree is and how it is used to control the movements in your structure and how to control it. You can use fold tree to also mention how different chains in a structure are supposed to move with respect to each other. For example, if your chain C is between chain A and B in a complex and its movements should affect chain B, you can use a fold tree that places chain C downstream of chain B. The same principles appply for any chains, including a ligand or a metal. So, you can control the behavior of their movements by applying different fold trees. If you are intrested in fold tree set up for a symmetruc pose, please check [here](https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/symmetry). Please note that a fold tree should contain NO CYCLEs. In other words, you cannot define a residue both as upstream and downstream of other set of residues. 
