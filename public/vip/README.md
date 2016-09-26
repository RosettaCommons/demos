# Rosetta VIP: Void Indentification and Packing

KEYWORDS: DESIGN STABILITY_IMPROVEMENT

This README was written in Feb. 2012, by Jim Havranek (havranek@genetics.wustl.edu).


This demo illustrates a protocol to identify candidate mutations for stabilizing proteins that have sub-optimally packed cores.

It has been published in the paper "Automated selection of stabilizing mutations in designed and natural proteins" by B. Borgo J.J. Havranek (2012), Proc. Natl. Acad. Sci. USA 109(5) pp 1494-99.

The source code is found at:  `(rosetta_source)/src/apps/public/vip.cc`

To generate the output in the example output directory, the app can be run as follows:

(where `ROSETTA3`=path-to-Rosetta/main/source)

```bash```
$> cp rosetta_inputs/* .
$> $ROSETTA3/bin/vip.default.linuxgccrelease @flags 
```

In the above, ".default.linuxgccrelease" may need to be replaced with your build, operating system, and compiler.  For example, for a static, debug-mode build on the Macintosh operating system with the clang compiler, you would use, ".static.macosclangdebug".

options are:

```
-cp:ncycles (Size)
	This will run the iterative protocol (find point mutations, relax, output best relaxed pose) a fixed number of times. If you don't use this option it will continue to run until it no longer finds favorable mutations. The latter option can take a while for large proteins.

-cp:max_failures (Size)
	This allows you to try each iteration multiple times before deciding it is pointless.  The reason for this option is that, due to RosettaHoles, the process is stochastic.  Whether a good mutation can be found can depend on the pseudorandom numbers pulled during the void identification steps.  The highest this should be set is around five.  The default is 1.

-cp:cutoff (real)
	This is the cutoff for choosing mutatable residues (ie distance from cavity ball to a non-bb, non-surface atom on the residue). The smallest cutoff you can use is best since that will mutated the smallest # of residues that line the cavity.

-sasa_calculator_probe_radius (real) Increasing this will likely give you surface clefts in addition to buried voids.
```

further options for the fast relax mover, sasa metric options, etc.:

```
-cp:pack_sfxn (score function) (e.g. -cp:pack_sfxn score12_full)
	Allows you to use a different score for the point mutant trials

-cp:relax_sfxn (score function) to use a different score for the relax stage

-cp:skip_relax
This causes the protocol to skip the relax step, which is a quick fix until I can exchange the relax mover for something more efficient. I'd only caution that scoring of the fixed bb point mutations isn't quantitative wrt to magnitude of the ddE. Rather, it seems to be very accurate in terms of favorable vs. unfavorable.

-cp:relax_mover
This controls the mover used for the relaxation step after inserting mutations.  The default is "relax", which performs the fast_relax protocol.  As available is "classic_relax", which can be a little slower.

-cp:local_relax
This restricts the scope of the relaxation protocol for possible mutations to 'neighbor' residues, which are defined as those residues within a X Ang. CB-CB distance.

-cp:print_reports
This generates a file (reports.txt ) that tells you which mutations were identified by the simple geometric screen, and which mutations still look good after relaxation.

-cp:print_intermediate_pdbs
This option outputs a fullatom pdb for each accepted mutation, named "vip_iter_1.pdb", with the number incremented for each pass through.  The pdb for the final iteration will be the same as the final pdb.

-cp:exclude_file
This allows you to specify a file that contains positions that should not be allowed to mutate.  The format is one position per line, with pdb number and chain separated by a space ("128 A").
```

a couple of other notes:

- It is preferable to let the application find as many mutations as possible.  This is accomplished by _not_ using the -cp:ncycles option.  If the application takes too long to run, the -cp:ncycles option can be used to limit the run time to an acceptable amount.

- Extra rotamers are automatically included in the point mutant trials, so if you use -ex1 -ex2 etc flags, these will be applied in the relax step which slows things down quite a bit.

- Cavity finding is done via RosettaHoles1, which is stochastic, so:  1) if it finds no cavities, run it again and it probably will and 2) separate runs can result in a different sequence of mutations.

