* *Description*

For this protocol, we will be using the [GenericMonteCarlo](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/GenericMonteCarloMover) 
mover, which is another type of "meta mover" that performs iterative trials of a given Mover and evaluates whether that move was *good* or *bad* based on a user-defined scoring metric.
If a move is *good*, the GenericMonteCarlo mover accepts or rejects that move based on the evaluation of a Boltzmann criterion.

Before we can use the GenericMonteCarloMover, we need to create a Mover that we can pass the GenericMonteCarloMover. Let's create a ParsedProtocol Mover that
first applies the SmallMover and then applies the PackRotamersMover. 

The SmallMover will make small-move style backbone torsion movements that do not minimize downstream propagation. 
The PackRotamersMover will alter the residue identities and/or sidechain rotamers of specific residues.

We would like to focus the SmallMover and PackRotamersMovers on a small subset of residues in the protein structure. To do so, let's add ResidueSelectors that will create four subsets of residues:
 * All residues in chain A (we will call this the "pdz" domain)
 * All residues in chain B (we will call this the "peptide")
 * All residues in the neighborhood of chain B (we will call this the "interface" of the peptide and domain)
 * All residues NOT in the interface region (we will call this "not interface")

In the RESIDUE_SELECTORS section of the script, add
```xml
  <RESIDUE_SELECTORS>
    <Chain name=pdz chains=A />  
    <Chain name=peptide chains=B />
    <Neighborhood name=interface selector=peptide distance=6 />
    <Not name=not_interface selector=interface/>  
  </RESIDUE_SELECTORS>
```

Now that we have grouped residues using the ResidueSelectors, we can turn the design off for the residues in the PDZ domain (chain A), and then turn off repacking as well for the residues in the PDZ domain that are not
in the interface between the PDZ domain and the peptide (chain B). In the TASKOPERATIONS section of the script, we can add

```xml
<TASKOPERATIONS>
  <InitializeFromCommandline name=init/>
  
  <OperateOnResidueSubset name=repackonly_pdz selector=pdz >
    <RestrictToRepackingRLT/>
  </OperateOnResidueSubset>
  
  <OperateOnResidueSubset name=fix_notinterface selector=not_interface >
    <PreventRepackingRLT />
  </OperateOnResidueSubset>

</TASKOPERATIONS>
```

The first task operation, InitializeFromCommandline, allows us to use the options `-ex1 -ex1 -use_input_sc` to influence the number of rotamers available to the Packer.

The second and third task operations both use the OperateOnResidueSubset task operation to apply Residue Level TaskOperations (RLTs) to the ResidueSelector named in the `selector` field.
The second task operation turns the design function of the Packer OFF for the residues selected in the `pdz` ResidueSelector.
The third task operation turns the design AND repacking function OFF for the residues selected in the `not_interface` ResidueSelector.

The combination of these task operations will result in chain B of our input structure being designable, the interface residues on chain A being repackable, and the rest of chain A being fixed.

Now, let's pass these ResidueSelectors and TaskOperations to the two Movers we set out to use. In the MOVERS section of the script, we are going to add
```xml
<MOVERS>

  <Small name=small_mover scorefxn=myscore temperature=0.5 nmoves=1 angle_max=6.0 preserve_detailed_balance=0 residue_selector=interface />
  <PackRotamersMover name=repack_mover scorefxn=myscore task_operations=init,repackonly_pdz,fix_notinterface />

</MOVERS>
```

> *Note that the SmallMover accepts ResidueSelectors, and the PackRotamersMover accepts TaskOperations.*

We can only pass one mover to the GenericMonteCarlo mover, and so to combine the SmallMover and the PackRotamersMover, we will use a ParsedProtocol meta-mover to apply these movers one after the other:

```xml
<MOVERS>

  <Small name=small_mover scorefxn=myscore temperature=0.5 nmoves=1 angle_max=6.0 preserve_detailed_balance=0 residue_selector=interface />
  <PackRotamersMover name=repack_mover scorefxn=myscore task_operations=init,repackonly_pdz,fix_notinterface />
  
  <ParsedProtocol name=small_repack mode=sequence >
    <Add mover_name=small_mover />
    <Add mover_name=repack_mover />
  </ParsedProtocol>

</MOVERS>
```

Now we are ready to build the GenericMonteCarlo Mover. To the MOVERS section, we will add:

```xml
  <GenericMonteCarlo name=gmc_mover mover_name=small_repack scorefxn_name=myscore sample_type=low temperature=0.8 trials=10 drift=1 preapply=false recover_low=1 />
```

In order to apply the GenericMonteCarlo Mover to our input structure, we need to add this mover to the PROTOCOLS section of the script:

```xml
<PROTOCOLS>
  <Add mover=gmc_mover />
</PROTOCOLS>
```



