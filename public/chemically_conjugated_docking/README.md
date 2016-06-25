Chemically Conjugated Docking
=============================

KEYWORDS: LIGANDS DOCKING

Updated by (parisah@uw.edu) to enable automated demo testing

Included are demos for three applications:

Step 1. UBQ_E2_thioester:

This demo contains the starting structures for the original published use of UBQ_E2_thioester, and is a copy of UBQ_E2_thioester's integration test. To run the demo (in the inputs subfolder):

```bash
$> cp UBQ_E2_thioester/inputs/* .
$> $ROSETTA3/bin/UBQ_E2_thioester.default.linuxgccrelease @options1
```

`$ROSETTA3`=path-to-Rosetta/main/source

Edit the options file first to set paths as needed.

Please refer to UBQ_E2_thioester's documentation, online or at rosetta_source/doc/apps/public/scenarios/UBQ_conjugated.dox, and also the publication, Saha A, Kleiger G, Lewis S, Kuhlman B, Deshaies RJ. Essential role for ubiquitin-ubiquitin-conjugating enzyme interaction in ubiquitin discharge from Cdc34 to substrate. Molecular Cell. 2011 Apr 8;42(1):75-83.

Note that the provided outputs are from the integration test, which runs in ~30 s. You will need to run the code for much longer (both longer individual runs, and many trajectories) to get scientifically useful results. The options file and documentation provide details on how to do that.

Step 2. UBQ_Gp_CYD-CYD:

This demo contains the starting structures for the original use of the UBQ_Gp series of executables, and is a copy of their integration tests. To run the demos (in the inputs subfolder):

```bash
$> cp UBQ_Gp_series/inputs/* .
$> $ROSETTA3/bin/UBQ_Gp_CYD-CYD.default.linuxgccrelease @options2
```

 Step 3: UBQ_Gp_LYX-Cterm

```bash 
$> $ROSETTA3/bin/UBQ_Gp_LYX-Cterm.default.linuxgccrelease @options3
```
Please refer to the documentation, online or at rosetta_source/doc/apps/public/scenarios/UBQ_conjugated.dox, and also the publication, Baker R, Lewis SM, Wilkerson EM, Sasaki AT, Cantley LC, Kuhlman B, Dohlman HG, Campbell SL. Site-Specific Monoubiquitination Activates Ras by Impeding GTPase Activating Protein Function. Submitted.

Note that the provided outputs are from the integration test, which runs in ~30 s. You will need to run the code for much longer (both longer individual runs, and many trajectories) to get scientifically useful results. The options file and documentation provide details on how to do that.

All are closely related.
These demos are copies of their integration tests.
See the linked READMEs for further details.
