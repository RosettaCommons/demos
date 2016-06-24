β-Strand Homodimer Design
=========================

KEYWORDS: DESIGN INTERFACES STRUCTURE_PREDICTION

This outlines how to use the applications involved in finding exposed beta-strands and then designing a protein with an exposed beta strand to be a homodimer. 
Written by Ben Stranges (stranges at unc dot edu)

There are three applications associated with this demo:

    rosetta/rosetta_source/src/apps/public/scenarios/beta_strand_homodimer_design/homodimer_design.cc
    rosetta/rosetta_source/src/apps/public/scenarios/beta_strand_homodimer_design/homodimer_maker.cc
    rosetta/rosetta_source/src/apps/public/scenarios/beta_strand_homodimer_design/exposed_strand_finder.cc

It would also be a good idea to look at the doxygen for these apps at:

    rosetta/rosetta_source/doc/apps/public/scenarios/beta_strand_homodimer_design.dox

It explains how everything works. 

The idea is that you scan through a list of pdbs and find ones with exposed beta strands, you then make potential homodimers along these strands then design the residues at the interface to stabilize the homodimer. 

Running the Demo
----------------

1.  Find exposed beta-strands:
    ```
    $> $ROSETTA3/bin/exposed_strand_finder.default.linuxgccrelease -s 2a7b_mpm.pdb.gz @finder_options > exposed_strands.txt
    ```

    The contents of `finder_options` are:

    * general options  
        `-ignore_unrecognized_res true`  
        `-packing::pack_missing_sidechains`  
        `-out::nooutput`: this protocol manages its own output, prevent job distributor from helping  
        `-mute core basic protocols.jd2.PDBJobInputter`

    * app specific options  
        `-beta_length 5`: how long of an exposed strand do you look for  
        `-sat_allow 2`: how many satisfied bb atoms do you allow in this range

    * allow alignment of a found strand to some target protein  
        `-check_rmsd false`: setting this to false prevents the code from doing rmsd comparisons

    * `-native anti_model.pdb`

    * `-strand_span B 5 11`

    In reality you'll probably want to use `-l` instead of `-s` to pass a bigger list of pdbs to look for exposed strands in.
    This is just an example that will work.
    When this runs look at the output in exposed_strands.txt, the important line is this one:

        ExposedStrand: FILE:         2a7b_mpm   CHAIN:    A   START:   806   END:   812   H_BONDS:   0

    This means that there is an exposed strand in pdb 2a7b_mpm in chain A between residues 806 and 812 and there are 0 satisfied bb_bb H bonds in every other residue along this span. 
    This information is then used in the next step.

    This application has another mode that allows you to match a found exposed strand onto a beta-strand involved in the interaction with another protein. It is still in development and not really known to work. Use at your own risk. To activate it you need to set -check_rmsd to true and pass a structure with -native and use the -strand_span option.

2.  Make the potential homodimers:  
    There are two ways to do this.
    If you are only doing it for one structure is is easy just to use this command line:(where `$ROSETTA3`=path-to-Rosetta/main/source)

    ```
    $> $ROSETTA3/bin/homodimer_maker.default.linuxgccrelease -s 2a7b_mpm.pdb.gz  @maker_options > maker_tracers
    ```
    ```
    -run::chain A
    -sheet_start 806
    -sheet_stop 812
    -window_size 5 
    -ignore_unrecognized_res true
    -mute core protocols.moves.RigidBodyMover basic.io.database
    ```

    If you need to run a bunch of these from the output of the exposed strand finder I have provided a script that reads a file (runner_input).
    This file should be multiple lines instead of just the one here.
    It's structure is
    ```
    /path/to/pdb/1av3.pdb chainletter betastart betaend
    ```
    You will need to modify some of the paths in runner.sh. Then run this command:
    ```
    ./runner.sh runner_input
    ```

    This outputs a bunch of pdbs:
    ```
    2a7b_mpm_A806_anti_wind_1_step_-1.pdb
    2a7b_mpm_A808_parl_wind_2_step_1.pdb
    2a7b_mpm_A808_parl_wind_2_step_0.pdb
    2a7b_mpm_A808_parl_wind_2_step_-1.pdb
    2a7b_mpm_A806_parl_wind_1_step_1.pdb
    2a7b_mpm_A806_parl_wind_1_step_0.pdb
    2a7b_mpm_A806_parl_wind_1_step_-1.pdb
    2a7b_mpm_A806_anti_wind_1_step_1.pdb
    ```

    However for these purposes we are only interested in `2a7b_mpm_A806_anti_wind_1_step_1.pdb` so I removed the rest in the interest of saving space.

3.  Homodimer design:  
    The next step is to take the output from above and make the files you need for symmetry.
    To do this you will need to use the symmetry script as so:
    ```
    $> perl $ROSETTA3/src/apps/public/symmetry/make_symmdef_file.pl -m NCS -a A -i B -p 2a7b_mpm_A806_anti_wind_1_step_1.pdb > symmdef
    ```

    This makes a bunch of files:
    ```
    2a7b_mpm_A806_anti_wind_1_step_1_symm.pdb
    2a7b_mpm_A806_anti_wind_1_step_1_model_AB.pdb
    2a7b_mpm_A806_anti_wind_1_step_1_INPUT.pdb
    2a7b_mpm_A806_anti_wind_1_step_1.kin
    symmdef
    ```

    However, the only one you really need is `2a7b_mpm_A806_anti_wind_1_step_1_INPUT.pdb` and `symmdef` so I removed the rest.
    Now you are ready for the full design runs. I suggest using mpi compiled executables but the command line below is general. 
    Run this command:
    ```
    $> $ROSETTA3/bin/homodimer_design.default.linuxgccrelease -s 2a7b_mpm_A806_anti_wind_1_step_1_INPUT.pdb.gz  -symmetry:symmetry_definition symmdef @design_options
    ```
    See the comments in the design_options file for descriptions of what does what.

    This will output designed structures with the name: `2a7b_mpm_A806_anti_wind_1_step_1_INPUT_000x.pdb.gz` and a score file: `score.fasc`.
    From there it is up to you to to chose how you will determine which designs meet your needs.

    Then you should probably run the output through the InterfaceAnalyzer. See documentation for it here:

    rosetta_source/doc/apps/public/analysis/interface_analyzer.dox
