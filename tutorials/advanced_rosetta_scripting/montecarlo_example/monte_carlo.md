* *Iteratively test backbone and sidechain movements.*
* *Record the trajectory of a protocol in a multimodel PDB file.*

For this protocol, we will be using the [GenericMonteCarlo](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/GenericMonteCarloMover) 
mover, which is another type of "meta mover" that performs iterative trials of a given Mover and evaluates whether that move was *good* or *bad* based on a user-defined scoring metric.
If a move is *good*, the GenericMonteCarlo mover accepts or rejects that move based on the evaluation of a Boltzmann criterion.

Let's consider the example of re-designing and refining a peptide that is bound to a protein structure in order to find a peptide with possibly higher binding affinity. In the figure below, the peptide is depicted in magenta sticks, whereas the protein domain is in rainbow cartoon with green sidechain sticks.

![2drk](https://github.com/RosettaCommons/demos/blob/XRW2016_kmb/tutorials/advanced_rosetta_scripting/montecarlo_example/figures/2drk.png)

We will allow the Packer to choose different identities and rotamers for the current peptide in its bound state, while also allowing the backbone torsion angles of the residues at the protein-peptide interface be sampled by a small number of degrees over many iterations. 

Before we can use the GenericMonteCarloMover for this task, we need to create a Mover that we can pass the GenericMonteCarloMover. Let's create a [ParsedProtocol](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/ParsedProtocolMover) Mover, which is another meta-over, which will first apply a [SmallMover](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/SmallMover) and then apply the [PackRotamersMover](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/PackRotamersMover). 

  * The SmallMover will make small-move style backbone torsion movements that do not minimize downstream propagation. 
  * The PackRotamersMover will alter the residue identities and/or sidechain rotamers of specific residues.

However, we also do not want to apply both movers to all residues, as that would be computationally expensive as well as result in changes to the protein in regions not directly related to the protein-peptide interface. Therefore, let's focus the SmallMover and PackRotamersMovers on a small subset of residues in the protein structure. To do so, we'll first add [ResidueSelectors](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ResidueSelectors) that will create four subsets of residues:
 * All residues in chain A (we will call this the "pdz" domain)
 * All residues in chain B (we will call this the "peptide")
 * All residues in the neighborhood of chain B (we will call this the "interface" of the peptide and domain)
 * All residues NOT in the interface region (we will call this "not interface")

In the `<RESIDUE_SELECTORS>` section of a new script, we will add the following lines:
```xml
...
  <RESIDUE_SELECTORS>
    <Chain name=pdz chains=A />  
    <Chain name=peptide chains=B />
    <Neighborhood name=interface selector=peptide distance=6.0 />
    <Not name=not_interface selector=interface/>  
  </RESIDUE_SELECTORS>
...
```
 * The [`Chain`](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ResidueSelectors#residueselectors_conformation-independent-residue-selectors_chainselector) residue selectors select all residues of a given chain. 
 * The [`Neighborhood`](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ResidueSelectors#residueselectors_conformation-dependent-residue-selectors_neighborhoodresidueselector) selector will select all residues that are a certain distance away from the residues in the `selector` selection, and will include the target residues as well. 
 * The [`Not`](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ResidueSelectors#residueselectors_logical-residueselectors_notresidueselector) selector reverses the selection specified in the `selection` tag and selects all resides that **not** in that selector.

Now that we have grouped residues using the ResidueSelectors, we can turn the design function off for the residues in the PDZ domain (chain A), and then turn off repacking as well for the residues in the PDZ domain that are not
in the interface between the PDZ domain and the peptide (chain B). In the `<TASKOPERATIONS>` section of the script, we can add the lines:

```xml
...
<TASKOPERATIONS>
  <InitializeFromCommandline name=init/>
  <OperateOnResidueSubset name=repackonly_pdz selector=pdz >
    <RestrictToRepackingRLT/>
  </OperateOnResidueSubset>
  <OperateOnResidueSubset name=fix_notinterface selector=not_interface >
    <PreventRepackingRLT />
  </OperateOnResidueSubset>
</TASKOPERATIONS>
...
```

The first task operation, [InitializeFromCommandline](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/InitializeFromCommandlineOperation), allows us to use the `-ex1 -ex1 -use_input_sc` commandline flags to influence the number of rotamers available to the Packer.

The second and third task operations both use the [OperateOnResidueSubset](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/OperateOnResidueSubsetOperation) task operation to apply [Residue Level TaskOperations](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/Residue-Level-TaskOperations) (RLTs) to the ResidueSelector named in the `selector` field.
* The second task operation turns the design function of the Packer OFF for the residues selected in the `pdz` ResidueSelector.
* The third task operation turns the design AND repacking function OFF for the residues selected in the `not_interface` ResidueSelector.

The combination of all of these task operations will result in chain B of our input structure being designable, the interface residues on chain A being repackable, and the rest of chain A being fixed.

Now, let's pass the ResidueSelectors and TaskOperations to the two Movers we set out to use. In the `<MOVERS>` section of the script, we will ad the lines:

```xml
...
<MOVERS>
  <Small name=small_mover scorefxn=myscore temperature=0.5 nmoves=1 angle_max=6.0 preserve_detailed_balance=0 residue_selector=interface />
  <PackRotamersMover name=repack_mover scorefxn=myscore task_operations=init,repackonly_pdz,fix_notinterface />
</MOVERS>
...
```

> *Notice that the SmallMover accepts ResidueSelectors, and the PackRotamersMover accepts TaskOperations.*

We can only pass one mover to the GenericMonteCarlo mover, and so to combine the SmallMover and the PackRotamersMover, we will use a [ParsedProtocol](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/ParsedProtocolMover) meta-mover to apply these movers in sequential order:

```xml
...
<MOVERS>

  <Small name=small_mover scorefxn=myscore temperature=0.5 nmoves=1 angle_max=6.0 preserve_detailed_balance=0 residue_selector=interface />
  <PackRotamersMover name=repack_mover scorefxn=myscore task_operations=init,repackonly_pdz,fix_notinterface />
  
  <ParsedProtocol name=small_repack mode=sequence >
    <Add mover_name=small_mover />
    <Add mover_name=repack_mover />
  </ParsedProtocol>
...
</MOVERS>
```

Now we are ready to build the GenericMonteCarlo Mover. To the `<MOVERS>' section, we will add:

```xml
...
  <GenericMonteCarlo name=gmc_mover mover_name=small_repack scorefxn_name=myscore sample_type=low temperature=0.8 trials=10 drift=1 preapply=false recover_low=1 />
...
```

This line will apply the two Movers in the ParseProtocol Mover to the input pose once, evaluate whether the score has reduced (`sample_type=low`), and if so, will accept this move with a Bolztmann probability that has a temperature factor of 0.8. If this move is accepted, the move is kept, and this new structure becomes the starting structure for the next GenericMonteCarlo trial (`drift=1`). This will continue until the total number of trials have complete (`trials=10`). At the end of the GenericMonteCarlo Mover, we can take the last structure as the output structure or the lowest-energy structure as the output structure. Here, we have specified to output the lowest-energy structure (`recover_low=1`). 

Don't forget to add the GenericMonteCarlo mover to the `<PROTOCOLS>` section of the script so that is applied to the pose:

```xml
<PROTOCOLS>
  <Add mover=gmc_mover />
</PROTOCOLS>
```

To visualize the trajectory the input pose hase taken during the GenericMonteCarlo mover, we can also add the [PBDTrajectoryRecorder](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/PDBTrajectoryRecorderMover)  Mover to the ParsedProtocols block. This will take snapshots of the pose as the moves are applied to it (regardless if they are accepted or not) and output them to a multimodel PDB file. So to the `<MOVERS>` section, let's add one more mover and add it to the ParsedProtocols block:

```xml

<MOVERS>
...
		<PDBTrajectoryRecorder name=pdb_traj_recorder stride=1 filename=traj.pdb cumulate_jobs=0 cumulate_replicas=0 />
		<ParsedProtocol name=small_repack mode=sequence >
		 <Add mover_name=small_mover />
		 <Add mover_name=repack_mover />
		 <Add mover_name=pdb_traj_recorder />
		</ParsedProtocol>
...
</MOVERS>
```

The output structure could be similar to this:

![2drk_post](https://github.com/RosettaCommons/demos/blob/XRW2016_kmb/tutorials/advanced_rosetta_scripting/montecarlo_example/figures/2drk_post.png)

Notice the changes in the peptide sequence from gray (input) to magenta, and also the small backbone and rotamer changes to the protein domain.

