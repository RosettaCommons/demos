KEYWORDS: NONCANONICALS GENERAL
This demo contains the starting structures for the original use of the UBQ_Gp series of executables, and is a copy of their integration tests.  To run the demos (in the inputs subfolder):

    UBQ_Gp_CYD-CYD.linuxgccrelease @options
    UBQ_Gp_LYX-Cterm.linuxgccrelease @options

Each uses the same input files.  Read the options file first to set local paths as needed.

Please refer to the documentation, online or at rosetta_source/doc/apps/public/scenarios/UBQ_conjugated.dox, and also the publication, Baker R, Lewis SM, Wilkerson EM, Sasaki AT, Cantley LC, Kuhlman B, Dohlman HG, Campbell SL.  Site-Specific Monoubiquitination Activates Ras by Impeding GTPase Activating Protein Function.  Submitted.

Note that the provided outputs are from the integration test, which runs in ~30 s.  You will need to run the code for much longer (both longer individual runs, and many trajectories) to get scientifically useful results.  The options file and documentation provide details on how to do that.
