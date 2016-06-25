Design the Rac/Raf interface
============================

KEYWORDS: DESIGN INTERFACES

This demo will walk through the steps of designing a protein-protein interface.
The goal of this protocol is to predict mutations on Raf that will enable it to 
bind Rac.

Demo files
----------

    README.md
    starting_files/
      -- 1c1y.pdb.gz
      -- 2ov2.pdb.gz
    rosetta_inputs/
      -- design_script.xml (dock/design script)
      -- raf-rac.pdb (raf-rac starting interface)
    scripts/
      -- score_vs_rmsd.R (creates score vs rmsd plot)
    output_files/
      -- directory to store output files

Running the Demo
----------------

1. Create model of Raf-Rac interaction.
    * Uncompress pdb files:
      ```
      gunzip 1c1y.pdb.gz
      gunzip 2ov2.pdb.gz
      ```
    * Open 1c1y and 2ov2 with Pymol.
    * Superimpose chain A of 1c1y with chain A of 2ov2.  Pymol command:
      ```
      super 2ov2 and chain A, 1c1y and chain A
      ```
    * Select modeled Raf-Rac complex. Pymol command:
      ```
      select raf-rac, (1c1y and chain B)+(2ov2 and chain A)
      ```
    * Save Molecule raf-rac.  In pymol:
      ```
      File->Save Molecule->raf-rac.pdb
      ```

2. Run dock/design protocol using rosetta_scripts.
    * Change working directory to the output directory.
      ```
      $> cd output_files
      ```
    * Run design_script.xml using rosetta_scripts executable.
      ```
      $> $ROSETTA3/bin/rosetta_scripts.macosgccrelease -s ../rosetta_inputs/raf-rac.pdb -parser:protocol ../rosetta_inputs/design_script.xml -in:file:native ../rosetta_inputs/raf-rac.pdb -ex1 -ex2 -ignore_unrecognized_res -nstruct 1 -overwrite
      ```
    * Options:
      * `-ex1 -ex2`: expand rotamer library for chi1 and chi2 angles used in repacking/design
      * `-ignore_unrecognized_res`: ignores HETATM lines in input PDB file
      * `-nstruct 1`: specifies how many times the protocol is run (one decoy is output for each run)
      * `-overwrite`: overwrites decoys from previous runs

    Each run should take about 90 seconds to complete.
    For production runs, nstruct should be set to 1000 or greater.
    This protocol returns decoys named `raf-rac_####.pdb` and a file named 
    `score.sc` that contains scores for each decoy.

3. Postprocessing output from dock/design run.
    * To see the sequence changes in chain B due to design grep out chain B: 
      grep " B " raf-rac.pdb > chainBin.pdb then use 
            rosetta_tools/protein_tools/scripts/SequenceProfile.py
      to diff the sequences

    * Run R script to generate Interface Score vs. Interface RMSD plot. R is 
      easy to install from [[here|http://cran.r-project.org/bin/macosx/]] as a 
      pkg that self installs (do "which R" to check install first). The above 
      link is fro MacOSX, but R is available for many other operating systems 
      as well.

            CMD BATCH scripts/score_vs_rmsd.R

      Output: `score_vs_rmsd.pdf`
