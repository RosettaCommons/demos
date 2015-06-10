AbInitio Structure Prediction Using Chemical-Shift Generated Fragments
======================================================================

Written by Lei Shi.
Ray Wang drafted the previous version.

---

We will use the chemical shifts to improve the fragments from which Rosetta builds up structures, and the NOEs to guide the Rosetta calculations towards the native structure.

Please see references at:
* rosetta abinitio: Bradley, P et al Science 2005
* chemical shift fragments: Shen Y et al. PNAS 2008;105:4685-4690
* chemical shift+NOE+RDC: Raman S, et al Science 2010

In this demo, we will use PDB 2JY7, which is a small protein (for demo purpose) and has experimental data deposited. Several scripts are provided in the scripts folder for formatting purposes:

    bmrb2talos.com
    cst_map_toCB.py
    upl2mini.csh
    scores.score.cfg

If you are from David Baker lab, there are scripts available to make setup easier without going through public servers. The following instructions should work just fine without having direct access to any Baker lab cluster.

Running the Demo
----------------
1. Create following folders:
    ```
    mkdir starting_inputs
    mkdir rosetta_inputs
    mkdir rosetta_inputs/talos_output
    mkdir rosetta_inputs/pick_cs_fragments
    ```

2. Download protein fasta and experimental data:  
Download fasta from http://www.pdb.org/pdb/explore/explore.do?structureId=2JY7
    ```
    wget http://www.pdb.org/pdb/files/fasta.txt?structureIdList=2JY7 -O starting_inputs/t000_.fasta
    ```
Download chemical shift data from http://www.bmrb.wisc.edu/data_library/summary/index.php?bmrbId=15591
    ```
    wget http://rest.bmrb.wisc.edu/bmrb/NMR-STAR2/15591 -O starting_inputs/raw.cs.bmrb
    ```

3. Format data for Rosetta use:  
Formmatting chemical shift data for TALOS:
    ```
    scripts/bmrb2talos.com starting_inputs/raw.cs.bmrb > rosetta_inputs/cs.talos
    ```

4. Generating talos predictions using http://spin.niddk.nih.gov/bax/nmrserver/talosn/ using `rosetta_inputs/cs.talos`.
Save/copy `pred.tab` and `predSS.tab` to `rosetta_inputs/talos_output`

5. Generate fragment/profile from RobettaServer http://www.robetta.org/fragmentqueue.jsp using `starting_inputs/t000_.fasta`.
Save/copy `t000_.checkpoint` to `rosetta_inputs/`

6. Pick fragments using secondary structure profile and chemical shift data:
    ```
    Rosetta/main/source/bin/fragment_picker -database Rosetta/main/database/ -in::file::vall Rosetta//tools/fragment_tools/vall.apr24.2008.extended.gz -frags::n_frags 200 -frags::frag_sizes 3 9 -frags::sigmoid_cs_A 2 -frags::sigmoid_cs_B 4 -out::file::frag_prefix rosetta_inputs/pick_cs_fragments/frags.score -frags::describe_fragments rosetta_inputs/pick_cs_fragments/frags.fsc.score -frags::scoring::config scripts/scores.score.cfg -in:file:fasta starting_inputs/t000_.fasta -in:file:checkpoint rosetta_inputs/t000_.checkpoint -in:file:talos_cs rosetta_inputs/cs.talos -frags::ss_pred rosetta_inputs/talos_output/predSS.tab talos -in::file::talos_phi_psi rosetta_inputs/talos_output/pred.tab
    ```

7. Run Rosetta with the fragments chemical shift fragments:
    ```
    Rosetta/main/source/bin/AbinitioRelax -database Rosetta/main/database/ -in:file:fasta starting_inputs/t000_.fasta -file:frag3 rosetta_inputs/pick_cs_fragments/frags.score.200.3mers -file:frag9 rosetta_inputs/pick_cs_fragments/frags.score.200.9mers -nstruct 1 -abinitio::increase_cycles 10 -abinitio::relax -score::weights score13_env_hb -abinitio::rg_reweight 0.5 -abinitio::rsd_wt_helix 0.5 -abinitio::rsd_wt_loop 0.5 -disable_co_filter true -out:file:silent csrosetta.out
    ```
Change nstruct to generate desired number of models. Larger is better depending on your available computer time, etc.

Processing the output
---------------------
1. Extract the low energy models:
    ```
    grep SCORE csrosetta.out | sort –nk2 | head
    ```
The second column contains the energies of the lowest energy 10 models.
Select as the cutoff the energy on the last line.

2. This script:
    ```
    cull_silent.pl csrosetta.out “score < cutoff”
    ```
will produce `csrosetta.select.silent` which contains the lowest models below cutoff.
  
3. Extract pdbs from selected silent file
    ```
    Rosetta/main/source/bin/extract_pdbs -database Rosetta/main/database/ -in::file::silent csrosetta.select.silent
    ```

4. Check convergence by superimposing the ten low energy models in pymol or your favorite molecular graphics

5. Check convergence by clustering the lowest energy models (see clustering demo for instructions)
