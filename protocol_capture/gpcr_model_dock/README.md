GPCR Model Dock Protocol Capture
================================
KEYWORDS: DOCKING GENERAL
Written in Jan 2013 by Elizabeth Nguyen (e dot dong dot nguyen at vanderbilt dot edu)

---

This demo contains all the files necessary to replicate the results from the 
paper:

Challenges of recovering native conformations of ligands docked into 
comparative models of G-protein coupled receptors. Nguyen,E.D.*, Norn, C.*, 
Frimurer, T., and Meiler,J. (2012), PLOS One. [[Supplementary 
Material|/protocol_capture/gpcr_model_dock/SupplementaryMaterial_ProtocolCapture.docx]].

Commands and files are provided to be able to replicate the following steps 
(from Figure 1):

1. Structural alignment of GPCR templates
2. Sequence alignment of the target GPCR to templates
3. Thread target sequence onto template backbone coordinates
4. Rebuild missing density
5. Rebuild ECL 1,2 and 3
6. Evaluate comparative models by clustering by full-receptor RMSD and 
   knowledge-based pocket residue filter
7. Generate ligand conformations in MOE
8. Dock ligand into comparative models
9. Analyze results by clustering binding modes by ligand RMSD


Files included in this protocol capture
---------------------------------------

Input files (rosetta_inputs):

    1u19A_clean.pdb 
    2vt4A_clean.pdb 
    2rh1A_clean.pdb 
    3emlA_clean.pdb 
    3oduA_clean.pdb 
    3pblA_clean.pdb 
    3rzeA_clean.pdb 
    3v2wA_clean.pdb 
    3uonA_clean.pdb 
    4dajA_clean.pdb 
    4dklA_clean.pdb 
    4djhA_clean.pdb 
    4ea3A_clean.pdb 
    4ej4A_clean.pdb
    all_gpcrs.fasta 
    1u19A.fasta
    1u19A.aln
    2rh1A_clean.pdb
    1u19A_on_2rh1A.pdb
    1u19A.jufo_ss
    1u19A.psipred_ss2
    1u19A.span
    1u19A.disulfide
    aa1u19A03_05.200_v1_3
    aa1u19A09_05.200_v1_3
    relax.options
    1u19A_on_2rh1A_relax.pdb
    ccd_initial.options
    1u19A_on_2rh1A.loops
    1u19A_on_2rh1A_initial.pdb
    1u19A_rmsd01.pdb
    1u19A.sdf
    1u19A.params
    1u19A_confs.pdb 
    1u19A_cluster01_01.pdb
    1u19A_cluster01_01_ligand.pdb
    dock.options
    dock.xml 
    1u19A_ligand.cluster.mat

Output files (example_outputs):

    1u19A_10percent_RMSD.txt
    cluster3_1u19A.Centers
    cluster3_1u19A.Rows 
    1u19A_cluster01_01_ligand_011u19A_cluster01_01_ligand_0001.pdb
    all.sdf
    cluster3_1u19A_ligand.Centers
    cluster3_1u19A_ligand.Rows

Scripts and applications not included in the Rosetta 3.4 release (scripts):

    evaluate_score_vs_pocket_rmsd
    jufo9d_span.pl
    rmsd.tcsh

Command lines to run this protocol capture
------------------------------------------

* Prepare GPCR crystal structures from the Protein Data Bank.

        rosetta_tools/protein_tools/scripts/clean_pdb.py 2RH1 A > 2rh1A_clean.pdb

* Perform a structural alignment of GPCRs using crystal structures from the Protein Data Bank. 

        mustang -p . -i 1u19A_clean.pdb 2vt4A_clean.pdb 2rh1A_clean.pdb 3emlA_clean.pdb 3oduA_clean.pdb 3pblA_clean.pdb 3rzeA_clean.pdb 3v2wA_clean.pdb 3uonA_clean.pdb 4dajA_clean.pdb 4dklA_clean.pdb 4djhA_clean.pdb 4ea3A_clean.pdb 4ej4A_clean.pdb -o all_gpcrs -F fasta -D 2.5 -s ON

* Sequence alignment of the target GPCR to templates

  Input target sequence 1u19A.fasta and profile alignment all_gpcrs.fasta to 
  http://mobyle.pasteur.fr/cgi-bin/portal.py#forms::clustalO-profile. We used 
  the default settings for this protocol.

* Thread target sequence onto template backbone coordinates

    rosetta_tools/protein_tools/scripts/thread_pdb_from_alignment.py --template=2rh1A_clean --target=1u19A --chain=A --align_format=clustal 1u19A.aln 2rh1A_clean.pdb 1u19A_on_2rh1A.pdb

* Generate secondary structure prediction, constraint file and fragments for bRh. 

  * Secondary structure- Jufo9D: http://meilerlab.org/index.php/servers/show?s_id=5

  * Secondary structure- PSIPRED: http://bioinf.cs.ucl.ac.uk/psipred/

  * Transmembrane span prediction based on Jufo9D:

            perl scripts/jufo9d_span.pl 1u19A.jufo9d > 1u19A.span 

  * Disulfide bond constraint file: Create file that lists residue number of 
    cysteine residues predicted to disulfide bond according to the alignment 
    with the template.

  * Fragment files: http://www.robetta.org

* Rebuilt missing density

        rosetta_source/bin/loopmodel.linuxgccrelease @ccd_initial.options -database rosetta_database 

* Rebuilt ECL 1,2 and 3 with CCD

        rosetta_source/bin/loopmodel.linuxgccrelease @ccd.options -database rosetta_database 

* Rebuilt ECL 1,2 and 3 with KIC

        rosetta_source/bin/loopmodel.linuxgccrelease @kic.options -database rosetta_database

* Analyze results by clustering top ten percent of comparative models by full receptor RMSD.

        bcl.exe PDBCompare -quality RMSD -atoms CA -pdb_list 1u19A_models.ls -aaclass AACaCb -prefix 1u19A_10percent_
        bcl.exe Cluster -distance_input_file 1u19A_10percent_RMSD.txt -input_format TableLowerTriangle -output_format Rows Centers -output_file cluster3_1u19A -linkage Average -remove_internally_similar_nodes 3

* Analyze results by filtering comparative models with a knowledge-based filter.

        scripts/evaluate_score_vs_pocket_rmsd/01_make_distances.csh
        scripts/evaluate_score_vs_pocket_rmsd/02_filter_models.py

* Create ligand conformations in MOE.

  See MOE operating guide. LowModeMD with the MMFFx94 force field and 
  Generalized Born solvation model was used to generate conformations within 
  the specified energy cutoff. The ligand conformations were then saved as an 
  .sdf file for conversion to .pdb and .params files for Rosetta.

        rosetta_source/src/python/apps/public/molfile_to_params.py -n 1u19A -p 1u19A 1u19A.sdf 
 
* Dock ligand into comparative models.

        rosetta_source/bin/rosettascripts.linuxgccrelease @dock.options -database rosetta_database

* Filter binding modes by energy, clustering and experimental restraints

        /scripts/rmsd.tcsh *.pdb
        bcl.exe ScoreSmallMolecule all.sdf output.sdf -comparison RMSD
        bcl.exe Cluster -distance_input_file 1u19A_ligand.cluster.mat -input_format TableLowerTriangle -output_format Rows Centers -output_file cluster3_1u19A_ligand -linkage Average -remove_internally_similar_nodes 3
