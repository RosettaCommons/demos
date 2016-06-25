AbInitio Structure Prediction Using Chemical-Shift Generated Fragments and NOE Distance Restraints
==================================================================================================

KEYWORDS: STRUCTURE_PREDICTION EXPERIMENTAL_DATA

=======
Written by Lei Shi.
Nikolas Sgourakis drafted the previous version.

---

We will use the chemical shifts to improve the fragments from which Rosetta builds up structures, and the NOEs to guide the Rosetta calculations towards the native structure.

Please see references at:
* rosetta abinitio: Bradley, P et al Science 2005
* chemical shift fragments: Shen Y et al. PNAS 2008;105:4685-4690
* chemical shift+NOE+RDC: Raman S, et al Science 2010

These Rosetta calculation steps are also described separately:
* Sgourakis NG et al JACS,2011,133(16):6288-98:

In this demo, we will use PDB 2JY7, which is a small protein (for demo purpose) and has experimental data deposited. Several scripts are provided in the scripts folder for formatting purposes:

	bmrb2talos.com
	cst_map_toCB.py
	upl2mini.csh
	scores.score.cfg

If you are from David Baker lab, there are scripts available to make setup easier without going through public servers. The following instructions should work just fine without having direct access to any Baker lab cluster.

Running the demo
----------------
1. Create following folders:  
    ```
    mkdir starting_inputs
    mkdir rosetta_inputs
    mkdir rosetta_inputs/talos_output
    mkdir rosetta_inputs/pick_cs_fragments
    ```

2. Download protein fasta and experimental data
Download fasta from http://www.pdb.org/pdb/explore/explore.do?structureId=2JY7  
    ```
    wget http://www.pdb.org/pdb/files/fasta.txt?structureIdList=2JY7 -O starting_inputs/t000_.fasta
    ```
Download chemical shift data from http://www.bmrb.wisc.edu/data_library/summary/index.php?bmrbId=15591  
    ```
    wget http://rest.bmrb.wisc.edu/bmrb/NMR-STAR2/15591 -O starting_inputs/raw.cs.bmrb
    ```
Download NOE data from http://restraintsgrid.bmrb.wisc.edu/NRG/MRGridServlet?pdb_id=2JY7&show_blocks=true&min_items=0:
    ```
    wget http://restraintsgrid.bmrb.wisc.edu/NRG/MRGridServlet?db_username=wattos1&format=ambi&mrblock_id=434910&pdb_id=2jy7&program=DYANA%2FDIANA&request_type=block&subtype=general+distance&type=distance
    echo "save file as starting_inputs/NOE_data.upl"
    ```

3. Format data for Rosetta use  
Formatting NOE: (Note only residues separated by more than 3 are kept in constraint)
The script `scripts/upl2mini.csh` only works with cyana format NOE:
    ```
    $> scripts/upl2mini.csh starting_inputs/NOE_data.upl > rosetta_inputs/NOE.cst
    $> scripts/cst_map_toCB.py rosetta_inputs/NOE.cst > rosetta_inputs/NOE.centroid.cst
    ```
Formmatting chemical shift data for TALOS. The required bmrb2talos script is part of the Rosetta toolbox, downloadable from [The CS Rosetta](http://csrosetta.chemistry.ucsc.edu/downloads/toolbox) web site. The output file is provided here.
    ```
    $ scripts/bmrb2talos.com starting_inputs/raw.cs.bmrb > rosetta_inputs/cs.talos
    ```

4. Generating talos predictions using http://spin.niddk.nih.gov/bax/nmrserver/talosn/ using rosetta_inputs/cs.talos
Save/copy pred.tab and predSS.tab to rosetta_inputs/talos_output

5. Generate fragment/profile from RobettaServer http://www.robetta.org/fragmentqueue.jsp using starting_inputs/t000_.fasta
Save/copy t000_.checkpoint to rosetta_inputs/

6. Pick fragments using secondary structure profile and chemical shift data: (where `$ROSETTA3`=path-to-Rosetta/main/source and `$ROSETTA3_DB`=path-to-database)
```
$> $ROSETTA3/bin/fragment_picker.default.linuxgccrelease -in::file::vall $ROSETTA3_DB/sampling/small.vall.gz -frags::n_frags 200 -frags::frag_sizes 3 9 -frags::sigmoid_cs_A 2 -frags::sigmoid_cs_B 4 -out::file::frag_prefix rosetta_inputs/pick_cs_fragments/frags.score -frags::describe_fragments rosetta_inputs/pick_cs_fragments/frags.fsc.score -frags::scoring::config scripts/scores.score.cfg -in:file:fasta starting_inputs/t000_.fasta -in:file:checkpoint rosetta_inputs/t000_.checkpoint -in:file:talos_cs rosetta_inputs/cs.talos -frags::ss_pred rosetta_inputs/talos_output/predSS.tab talos -in::file::talos_phi_psi rosetta_inputs/talos_output/pred.tab
```
**IMPORTANT**
The *small.vall.gz* used here for fragment picking is only used to speed up the demo. You have to change this to the vall database on your system!

7. Run Rosetta with the fragments made above and use NOEs to guide search
```
$> ROSETTA3/bin/minirosetta.default.linuxgccrelease -cst_fa_file rosetta_inputs/NOE.cst -cst_file rosetta_inputs/NOE.centroid.cst -abinitio:stage1_patch scripts/patch_atom_pair_constraint -abinitio:stage2_patch scripts/patch_atom_pair_constraint -abinitio:stage3a_patch scripts/patch_atom_pair_constraint -abinitio:stage3b_patch scripts/patch_atom_pair_constraint -abinitio:stage4_patch scripts/patch_atom_pair_constraint -score:patch scripts/patch_atom_pair_constraint -in:file:fasta starting_inputs/t000_.fasta -file:frag3 rosetta_inputs/pick_cs_fragments/frags.score.200.3mers -file:frag9 rosetta_inputs/pick_cs_fragments/frags.score.200.9mers -nstruct 1 -out:file:silent csrosetta_noe.out -run:protocol abrelax -abinitio::increase_cycles 0.1 -overwrite -abinitio::relax
```
**IMPORTANT**  
* Change nstruct to generate the desired number of models. Larger is better depending on your available computer time, etc.
* Change -abinitio::increase_cycles 0.1 to 10! (0.5 is chosen for tesing purposes only)


You can/should adjust the weights of NOE constraints in `scripts/patch_atom_pair_constraint`.
You should also change nstruct to generate desired number of models.
Larger is better depending on your available computer time, etc.
Note that the demo in abinitio_w_chemicalshift_only, you can add flags such as:
```
    -abinitio::rg_reweight 0.5
    -abinitio::rsd_wt_helix 0.5
    -abinitio::rsd_wt_loop 0.5
    -disable_co_filter true
    -abinitio::increase_cycles 10
```
to help sampling in centroid stage.
They are not used probably NOE constraint helps guided the search.

Processing the output
---------------------
1. Extract the low energy models:
    ```
    grep SCORE csrosetta_noe.out | sort –nk2 | head
    ```
The second column contains the energies of the lowest energy 10 models.
Select as the cutoff the energy on the last line.
You should also use NOE constraint energy as a criteria to select structures.
Example is only provided for total score.

2.  This command
    ```
    cull_silent.pl csrosetta_noe.out “score < cutoff”
    ```
will produce csrosetta.select.silent which contains the lowest energy 10 models.

3. Extract pdbs from selected silent file
    ```
    Rosetta/main/source/bin/extract_pdbs -in::file::silent csrosetta.select.silent
    ```

4. Check convergence by superimposing the ten low energy models in pymol or your favorite molecular graphics.

5. Check convergence by clustering the lowest energy models (see clustering demo for instructions).

6. To see how NOE constraints are satisfied by a model: (lowscore_1.pdb is just an example and you can replace it with your best output)
    ```
    $> $ROSETTA3/bin/r_cst_tool.linuxgccrelease -in:file:s rosetta_inputs/lowscore_1.pdb -cst_file rosetta_inputs/NOE.cst
    ```
`r_cst_tool` is a pilot program by Oliver Lange in `Rosetta/main/source/src/apps/pilot/olli/`
