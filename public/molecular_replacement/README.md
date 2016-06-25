Molecular Replacement Demo
==========================

KEYWORDS: EXPERIMENTAL_DATA STRUCTURE_PREDICTION

Written October 26, 2010
Modified August 27, 2013
Modified June 24, 2016

---

This document briefly walks through the use of Rosetta to solve difficult
molecular replacement problems.  These tools assume that the user has access to
the Phenix suite of crystallographic software (in particular, phaser and the
mapbuilding script mtz2map); however, all intermediate files are included so
that if the user does not, most of the demo may still be run.

The basic protocol is done in 5 steps; each step has a corresponding script in
the folder:

1. Using HHSearch, find potential homology to the target sequence.  Use a
   Rosetta "helper script" to prepare templates (and Rosetta inputs for
   subsequent computations).

2. Use PHASER to search for placement of the trimmed templates within the unit
   cell.

3. Generate a map correspoding to each putative MR solution.

4. Using Rosetta, rebuild gaps and refine each template/orientation in Rosetta,
   constrained by the density of each solution.  After rescoring with PHASER,
   the best template/orientation should be clear (if the correct solution was
   among the starting models).

Step 1: prepare\_template\_for\_MR.sh
----------------------------------

This command-line illustrates the use of my script for preparing templates for
an initial phaser run.  Functionally, it's doing the same thing as the
crystallographic software 'Sculptor' but it doesn't remap the residues as
sculptor does (and makes it easier to run with different alignments).  The
script takes just one arguments: an HHR format alignment file.

Alignments generally come from HHsearch's web interface
(http://toolkit.tuebingen.mpg.de/hhpred).  After submitting the sequence
through their website, export the results to a .hhr file.  Results may be
trimmed so only alignments with a reasonable e-value and sequence coverage are
included.

The script parses the .hhr file, downloads each template PDB, and trims the PDB
to the aligned residues.  In addition, the script produces a 'rosetta-style'
alignment file; the format is briefly introduced below.  These alignment files
are used in Rosetta model-building.

    ## 1CRB_ 2qo4b
    # hhsearch
    scores_from_program: 0 1.00
    2 DFNGYWKMLSNENFEEYLRALDVNVALRKIANLLKPDKEIVQDGDHMIIRTLSTFRNYIMDFQVGKEFEEDLTGIDDRKCMTTVSWDGDKLQCVQKGEKEGRGWTQWIEGDELHLEMRAEGVTCKQVFKKV
    0 AFSGTWQVYAQENYEEFLRAISLPEEVIKLAKDVKPVTEIQQNGSDFTITSKTPGKTVTNSFTIGKEAEIT--TMDGKKLKCIVKLDGGKLVCRTD----RFSHIQEIKAGEMVETLTVGGTTMIRKSKKI
    --

The first line is '##' followed by a code for the target and one for the
template.  The second line identifies the source of the alignment; the third
just keep as it is.  The fourth line is the target sequence and the fifth is
the template ... the number is an 'offset', identifying where the sequence
starts.  However, the number doesn't use the PDB resid but just counds residues
_starting at 0_.  The sixth line is '--'.

The results for this demo appear in the folder 'templates'.  For each
alignement in the starting .hhr file, 3 files are produced.

You can run the file either by running the provided .sh file or:

```
$> $ROSETTA3/src/apps/public/electron_density/prepare_template_for_MR.pl inputs/1crb.hhr
```
where `$ROSETTA3`=path-to-Rosetta/main/source

Steps 2 & 3: run\_phaser.sh and make\_maps.sh
-------------------------------------------

This command line shows the use of Phaser to generate initial molecular
replacement solutions.  For each template we run phaser to find potential
placements of each template in the unit cell.

*NOTES* these steps require havin [PHENIX](https://www.phenix-online.org/) installed.

The example scripts here only generate a single model from a single template,
but for a real-world case, one will often want to use many different templates
and may want to generate more than one possible solution using 'TOPFILES n'.
In general, though, we have found it is better to use fewer potential solutions
from more templates than many solutions from few templates.

Sometimes weak hits may be found by lowering the rotation function cutoff in
phaser by adding the line 'SELECT ROT FINAL PERCENT 0.65' (or even 0.5) to the
phaser script.  Increasing the packing function threshold (with PACK 10) may
also help in some cases.

Finally, for each template/orientation, we generate the 2mfo-dfc map for input
to Rosetta in the next step.

Steps 4A & 4B: run\_rosetta\_mr.sh
--------------------------------

The final step illustrate the use of rosetta's comparative modeling into
density.  After running the script and an initial phaser run, density maps are
generated from each phaser hit, and cm-into-density is done.  The flag
-MR::mode cm is used to run this mode.  This first application does not try to
rebuild gaps in the alignment, it just performs the threading and runs relax
into density.  Thus, the only inputs needed are: the target fasta file, the
rosetta-style ali file, and the template pdb.  Because there is no rebuilding,
not many models are needed to adequately cover conformational space, generally
10-20 is sufficient.

This script is the same as above, but also rebuilds gaps in the alignment.  The
main difference is that a non-zero value is given for
`-MR::max_gaplength_to_model`; additionally, some flags must be given that
describe how rosetta should rebuild gaps.

<<<<<<< HEAD
Several additional input files must be provided as well.  Rebuilding of gaps is
done by fragment insertion (as in Rosetta ab initio); thus two backbone
fragment files (3-mers and 9-mers) must be given.  The application for building
these is included with rosetta but requires a bunch of external
tools/databases.  The easiest way to generate fragments is to use the Robetta
server (http://robetta.bakerlab.org/fragmentsubmit.jsp).  The fragment files
should be built with the full-length sequence; rosetta handles remapping the
=======
Several additional input files must be provided as well.  Rebuilding of gaps is 
done by fragment insertion (as in Rosetta ab initio); thus two backbone 
fragment files (3-mers and 9-mers) must be given.  The application for building 
these is included with rosetta but requires a bunch of external 
tools/databases.  The easiest way to generate fragments is to use the Robetta 
server [[http://robetta.bakerlab.org/fragmentsubmit.jsp](http://robetta.bakerlab.org/fragmentsubmit.jsp).  The fragment files 
should be built with the full-length sequence; rosetta handles remapping the 
>>>>>>> cbbf46a34eabb7ee1743531183eb99e185343849
fragments if not all gaps are rebuilt.

A brief overview of flags is given below:

```
    -in::file::fasta inputs/1crb.fasta
    -in::file::alignment inputs/1crb_2qo4.ali
    -in::file::template_pdb inputs/2qo4.PHASER.1.pdb
```
The fasta, alignment and template PDBs.  See section 1 for the input file format if it needs to be hand-edited.

```
    -edensity:mapfile inputs/sculpt_2QO4_A.PHASER.1.map
    -edensity:mapreso 3.0
    -edensity:grid_spacing 1.5

```

This is how the density map and scorefunction parameters are given to Rosetta.  The input map (-edensity:mapfile) is CCP4 format.  The flags 'mapreso' defines the resolution of the calculated density; If the data is high-resolution it is often good to limit this to 2.5 or 3.  The grid spacing should be no more than 1/2 the map resolution; if -MR::fast is used (see below), then the grid_spacing flag may be omitted (since more finely sampled grids will have much less of a speed penalty).

```
    -MR::cen_dens_wt 4.0
    -MR::fa_dens_wt 1.0
```

This controls the weight on the experimental density data furing the two stages of refinement.  The second flag (fa_dens_wt) has the greatest impact on the final models.  If, after model generation and visual inspection, models don't seem to be fitting the density well (or overfitting to the density), this may be adusted accordingly.  If omitted, the values shown abve are the defaults that are used; generally, these defaults are sufficient for many cases.

   ` -MR::fast`   
        A special faster density scoring formulation is used.  Off by default, but it is recommended.

   ` -MR::max_gaplength_to_model 8`  
        Rosetta will close gaps up to this width; the larger this value is, the more sampling is required.  Values higher than 10 will often return incorrect loop conformations, although for very restrained segments, or largely helical segments, large insertions may be successfully modeled.

    -nstruct 20   
   The number of output structures.  Generally 10-20 is sufficient, unless a large 'max_gaplength_to_model' is given.

    -ignore_unrecognized_res   
   If the template contains nonstandard residues/ligands/waters, this tells Rosetta to ignore them.  This flag is recommended.

    -loops::frag_files inputs/aa1crb_09_05.200_v1_3.gz inputs/aa1crb_03_05.200_v1_3.gz none   
   (Optional) 	Fragment files from Robetta.  If omitted, MR-Rosetta will automatically generate fragments for the input structure; this may slightly reduce final model accuracy.  (The two separate command lines illustrate using and omitting this flag).

-

Since each model is independently generated, multiple processes may be used to
produce all the necessary models.  To manage the output, either each process
can be run from a separate directory, or '-out:prefix <prefix>' can be used to
keep jobs from overwriting each other's structures.  Rosetta workloads may also
be split using MPI; see the rosetta documentation for more details.

For a short test of these, you can run these commands for steps 4 and 5, respectively, using provided inputs. Please note that for Step 4 you need to generate fragments (by submitting your pdb to http://robetta.bakerlab.org/).
```
$> $ROSETTA/bin/mr_protocols.default.linuxclangrelease -in::file::fasta inputs/1crb.fasta -in::file::alignment templates/2qo4.ali -in::file::template_pdb phaser/2qo4_mr.PHASER.1.pdb -loops::frag_files inputs/frags.200.3mers inputs/frags.200.3mers none -edensity:mapreso 3.0 -edensity:grid_spacing 1.5 -edensity:mapfile phaser/2qo4_mr.PHASER.1_2mFo-DFc.ccp4 -MR::max_gaplength_to_model 8 -MR::fast -nstruct 1 -ignore_unrecognized_res -overwrite

```
$> $ROSETTA/bin/mr\_protocols.default.linuxclangrelease -in::file::fasta inputs/1crb.fasta -in::file::alignment templates/2qo4.ali -in::file::template\_pdb phaser/2qo4\_mr.PHASER.1.pdb -edensity:mapreso 3.0 -edensity:grid\_spacing 1.5 -edensity:mapfile phaser/2qo4\_mr.PHASER.1\_2mFo-DFc.ccp4 -MR::max\_gaplength\_to\_model 8 -MR::fast -nstruct 1 -ignore\_unrecognized\_res -overwrite
```

Alternatively, there is a compact output format, 'silent files' (see [[Controlling Input and Output in Rosetta]]) that can be used to dump structures to.  Simply add the flags '-out:file:silent <silent_filename> -out:file:silent_struct_type binary' and all structures from one process will be written to this compact file.  Then the rosetta program 'extract_pdb' can be used to extract:

    $ bin/extract_pdbs.default.linuxgccrelease -database $DB -in:file:silent <silent_filename> -silent_struct_type binary

