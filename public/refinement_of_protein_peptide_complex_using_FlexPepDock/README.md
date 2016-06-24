Refinement of peptide-protein complexes using FlexPepDock refinement
--------------------------------------------------------------------

KEYWORDS: STRUCTURE_PREDICTION PEPTIDES INTERFACES DOCKING

This demo illustrates how to run FlexPepDock refinement of a peptide-protein complex. The FlexPepDock Refinement protocol is designed to create high-resolution models of complexes between a flexible peptide and a globular protein, with side chains of binding motifs modeled at nearly atomic accuracy. It is intended for cases where an approximate, coarse model of the interaction is available.

Protocol Overview
-----------------
The input to the refinement protocol is a coarse model of the peptide-protein complex in PDB format. The protocol iteratively optimizes the peptide backbone and its rigid-body orientation relative to the receptor protein, including periodic on-the-fly side-chain optimization. The protocol is able to account for a considerable diversity of peptide conformations within a given binding site. However it is important to note that the refinement protocol can only refine models which are close to the correct solution both in terms of Cartesian and dihedral (phi/psi) distance. For the cases where no information regarding the peptide backbone is availavle, we recommend to use the FlexPepDock ab-initio protocol (See the demo 'abinitio_fold_and_dock_of_peptides_using_FlexPepDock').


Running the FlexPepDock refinement protocol:
--------------------------------------------
1. Create an initial complex structure: A coarse model of the peptide-protein interaction can be obtained from a low resolution peptide-protein docking protocol or can be built using a homologous structure. Here in this demo we have provided 1AWR.ex.pdb in which the peptide is in extended conformation. This will serve as the coarse model and our goal will be to refine it to obtain a near-native model. The native structure is 1AWR.pdb. Both 1AWR.ex.pdb and 1AWR.pdb are located in the input directory.

2. Prepack the input model: This step involves the packing of the side-chains in each monomer to remove internal clashes that are not related to inter-molecular interactions. The prepacking guarantees a uniform conformational background in non-interface regions, prior to refinement. The prepack_flags file contains the flags for running the prepacking job. The run_prepack script will run prepacking of the input structure 1AWR.ex.pdb located in the input directory.

Run the prepack as follows: (`$ROSETTA3`= path-to-Rosetta/main/source)
```
$> $ROSETTA3/bin/FlexPepDocking.default.linuxgccrelease @prepack_flags >output/ppk.log
$> mv input/1AWR.ex_0001.pdb input/1AWR.ex.ppk.pdb
```

The output will be a prepacked structure, 1AWR.ex.ppk.pdb located in the input directory; a scorefile named ppk.score.sc and a log file named prepack.log file located in the output directory. This prepacked structure will be used as the input for the refinement step.

3. Refine the prepacked model: This is the main part of the protocol. In this step, the peptide backbone and its rigid-body orientation are optimized relative to the receptor protein using the Monte-Carlo with Minimization approach, including periodic on-the-fly side-chain optimization. An optional low-resolution (centroid) pre-optimization will increase the sampling range and may improve performance further. The file refine_flags contains flags for running the refinement job. The run_refine script will run refinement of the prepacked structure generated in the prepacking step located in the input directory.

Run the refine using the following commands:

```
$> $ROSETTA3/bin/FlexPepDocking.default.linuxgccrelease  @refine_flags >refine.log
echo "FlexPepDock refinement done!"
```

The output will be a refined structure (1AWR.ex.ppk_0001.pdb) located in the output directory; a scorefile named refine.score.sc and a log file named refine.log file located in the output directory. This script has to be modified to run on a cluster during a production run. You may see an warning in the tracer saying "Inaccurate G! step= 3.8147e-06 Deriv= -0.0101137 Finite Diff= 0.00449653". No worries, that is a known Rosetta bug.

Further information
-------------------
A detailed documentation on FlexPepDock is available at https://www.rosettacommons.org/docs/latest/application_documentation/docking/flex-pep-dock
Raveh B, London N, Schueler-Furman O. (2010).Sub-angstrom modeling of complexes between flexible peptides and globular proteins. Proteins 78:2029-40.

