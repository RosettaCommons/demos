# Creating protocols with RosettaScripts
======================================

KEYWORDS: SCRIPTING_INTERFACES CORE_CONCEPTS

Tutorial by Rocco Moretti (rmorettiase@gmail.com) and Vikram K. Mulligan (vmullig@uw.edu).  Created on 21 June 2016 as part of the 2016 Documentation XRW.

Updated 29 May 2017 by Vikram K. Mulligan (vmullig@uw.edu) for new ref2015 scorefunction.

[[_TOC_]]

## Goals

At the end of this tutorial, you will understand:

- The RosettaScripts paradigm
- RosettaScripts syntax
- How to control final file output from RosettaScripts
- How to manipulate poses in RosettaScripts using *movers*
	- How to control movers that invoke the minimizer using *MoveMaps*
	- How to control movers that invoke the packer using *TaskOperations*
- How to select residues in RosettaScripts using *residue selectors*
- How to evaluate pose properties and control protocol flow in RosettaScripts using *filters*

## What is RosettaScripts?

Originally, the interface for Rosetta3 functionality was individual applications,
each made specifically for a particular use. One drawback of this approach was that
customization of protocols was difficult. If the protocol author properly anticipated
users' needs, then they may have put in options which allowed users to change the protocol.

However, these applications are typically limited in the extent to which they allow users
to modify the protocol. To allow for greater flexibility, RosettaScripts was created.
RosettaScripts allows users to create and modify protocols using an XML based syntax.  Broadly,
RosettaScripts is based around the paradigm of having a single structure (the *pose*) that
enters the protocol, a series of steps performed, each modifying the pose in some way (*movers*)
or evaluating some property of the pose (*filters*), and a single structure written out.  The
protocol can then be run repeatedly to generate large ensembles of output structures, or to
process large ensembles of input structures.  Even more broadly, RosettaScripts lets a user link
individual Rosetta modules together in a linear sequence.

This tutorial is intended to take you through the process of creating a new protocol
with RosettaScripts. It should also give you a good grounding in how you can modify
existing RosettaScripts protocols. Note that you can certainly run RosettaScripts without
modifying the XML - the most common use case of RosettaScripts is probably re-using
an XML produced by someone else.

## Your first RosettaScript

* *Run the simplest possible RosettaScript*

The simplest RosettaScript XML is one which does nothing.  You can obtain a skeleton XML file, which does nothing, in one of two ways.  You go to [the RosettaScripts documentation page](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/RosettaScripts) and find the skeleton XML file there, then copy and paste it into a new file (`nothing.xml`).  You can also generate a skeleton XML file by running the rosetta_scripts application without any parameters:

```bash
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease > nothing.xml
```

In the above, the ".default.linuxgccrelease" may need to be changed for your build, operating system, and compiler (*e.g.* ".static.macosclangrelease" for the static build using the clang compiler on the Macintosh operating system).  If you run the above, it will produce output similar to the following:

```
core.init: Rosetta version unknown:8c0acaa002dc2aeea2b97370cad3a5c9e1f6b2fe 2016-09-22 12:17:48 -0500 from git@github.com:RosettaCommons/main.git
core.init: command: /ssd1/morettr/Rosetta6/main/source/bin/rosetta_scripts.default.linuxgccrelease
core.init: 'RNG device' seed mode, using '/dev/urandom', seed=964549876 seed_offset=0 real_seed=964549876
core.init.random: RandomGenerator:init: Normal mode, seed=964549876 RG_type=mt19937
core.init: Resolved executable path: /ssd1/morettr/Rosetta6/main/source/build/src/release/linux/2.6/64/x86/gcc/5.2/default/rosetta_scripts.default.linuxgccrelease
core.init: Looking for database based on location of executable: /ssd1/morettr/Rosetta6/main/database/
core.init:
core.init: USEFUL TIP: Type -help to get the options for this Rosetta executable.
core.init:
apps.public.rosetta_scripts.rosetta_scripts: No XML file was specified with the "-parser:protocol <filename>" commandline option.  In order for RosettaScripts to do something, it must be provided with a script.
apps.public.rosetta_scripts.rosetta_scripts: The following is an empty (template) RosettaScripts XML file:

<ROSETTASCRIPTS>
	<SCOREFXNS>
	</SCOREFXNS>
	<RESIDUE_SELECTORS>
	</RESIDUE_SELECTORS>
	<TASKOPERATIONS>
	</TASKOPERATIONS>
	<FILTERS>
	</FILTERS>
	<MOVERS>
	</MOVERS>
	<APPLY_TO_POSE>
	</APPLY_TO_POSE>
	<PROTOCOLS>
	</PROTOCOLS>
	<OUTPUT />
</ROSETTASCRIPTS>

At any point in a script, you can include text from another file using <xi:include href="filename.xml" />.
apps.public.rosetta_scripts.rosetta_scripts: Variable substituion is possible from the commandline using the -"parser:script_vars varname=value" flag.  Any string of the pattern "%%varname%%" will be replaced with "value" in the script.
apps.public.rosetta_scripts.rosetta_scripts:
apps.public.rosetta_scripts.rosetta_scripts: The rosetta_scripts application will now exit.
```

This will be written to the standard output, which has been redirected to the file nothing.xml.
You can delete all lines preceding ```<ROSETTASCRIPTS>``` and following ```</ROSETTASCRIPTS>``` to obtain a minimal template.

Before running this script, let's edit it slightly to add comments:

```xml
<ROSETTASCRIPTS>
	<SCOREFXNS>
	</SCOREFXNS>

	This is a comment

	<RESIDUE_SELECTORS>
	</RESIDUE_SELECTORS>
	<TASKOPERATIONS>
		           So is this
	</TASKOPERATIONS>
	<FILTERS>
	</FILTERS>
	<MOVERS>

	  Anything not in angle brackets* is a comment.
	  This makes it easy to temporarily disable things by deleting just the first angle bracket.

	  MyMover name=mover1 option1="false" option2="23" /> Here is a mover that is commented out and ignored by RosettaScripts.  If I add back an angle bracket before "MyMover", it will be parsed.

	</MOVERS>
	<APPLY_TO_POSE>
	</APPLY_TO_POSE>
	<PROTOCOLS>
	</PROTOCOLS>
	<OUTPUT />

*(Angle brackets are the greater than/less than signs)

</ROSETTASCRIPTS>
```

The nothing.xml file is also provided in the inputs directory:

```bash
$> cp inputs/nothing.xml .
```

As you haven't further defined any protocol, this XML does nothing to the structure. As a test, let's just run a structure through RosettaScripts with this XML. RosettaScripts takes the standard input and output flags. In addition, the `-parser:protocol` option specifies which XML file to use.

```bash
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol nothing.xml
```

In the tracer output, Rosetta should print its interpretation of the XML input.

```xml
<ROSETTASCRIPTS>
	<SCOREFXNS/>
	<RESIDUE_SELECTORS/>
	<TASKOPERATIONS/>
	<FILTERS/>
	<MOVERS/>
	<APPLY_TO_POSE/>
	<PROTOCOLS/>
	<OUTPUT/>
</ROSETTASCRIPTS>
```

The first thing to notice is that the comments added to the XML (everything outside the angle brackets) is ignored.

Secondly, this demonstrates different ways of writing XML tags. XML tags are surrounded by angle brackets (greater/less than signs). A tag must be closed by a slash. Tags can be nested in other tags (like SCOREFXNS is nested within ROSETTASCRIPTS), in which case the outer tag must be closed by something like `</ROSETTASCRIPTS>`. If the tags are not nested, they can be closed by putting the slash at the end of the tag, like `<SCOREFXNS/>`.  The following two statements are perfectly equivalent:

```xml
<SCOREFXNS>
</SCOREFXNS>
```
```xml
<SCOREFXNS/>
```

In the above, the latter is more concise, though, at the expense of preventing anything from being enclosed within the SCOREFXNS block.

Additionally, whitespace is largely ignored in RosettaScripts.  The following three statements are perfectly synonymous:
```xml
<SCOREFXNS></SCOREFXNS>
```
```xml
<SCOREFXNS>      </SCOREFXNS>
```
```xml
<SCOREFXNS>
</SCOREFXNS>
```

Conventionally, tags are indented in proportion to their level of nesting, but this is for human readability, not for machine parsing; the rosetta\_scripts application disregards tabs entirely.  The one case in which whitespace matters is when setting options within a tag.  When a tag contains an option that accepts a comma-separated list, these must *not* have whitespace within them:

```xml
<PackRotamers name="pack1" task_operations="task1,task2,task3" /> #This is allowed
<PackRotamers name="pack2" task_operations="task2, task2, task3" /> #This will be misinterpreted
```

This brings up another RosettaScripts syntax convention: generally, we have blocks that define *types* of objects, and within these blocks, we define individual *instances* of objects of the type, giving each one a unique name.  For example, the ```<MOVERS> ... </MOVERS>``` block is the place to define movers.  Within this, we define specific instances of specific types of movers, and we set options for these movers, including a unique name by which each mover will be addressed at later points in the script.  For example:

```xml
	<MOVERS>  #In this section, movers are defined.
		  #The following is a particular mover of the "PackRotamers" type, which we give the
		  #unique name "pack1".  It takes, as an option, a list of previously-defined
		  #TaskOperation objects (a type of object that will be introduced later in this
		  #tutorial).  We assume that task1, task2, and task3 were defined and given these
		  #unique names prior to this point in the script.

		  <PackRotamers name="pack1" task_operations="task1,task2,task3" />

		  #From now on, we can refer to the mover defined above using the unique name "pack1".
	</MOVERS>
```

Looking at the output PDB, the output structure (1ubq\_0001.pdb) should be nearly identical to the input structure. The major difference should be the presence of hydrogens which were not in the input structure. This is *not* something that is specific to RosettaScripts - in general Rosetta will add missing hydrogens and repack sidechain atoms missing in the input PDB.

Additionally, you should see the standard Rosetta score table at the end of the PDB. By default, the structure will be rescored with the default Rosetta score function (ref2015, as of this writing). This can be controlled by the ```-score:weights``` command line option.

## Controlling RosettaScripts File Output

* *Score the output with a custom scorefunction*
* *Control output file format*

Before we explore the full power of RosettaScripts, let's make sure that we understand how to control the rosetta\_scripts application's output.  There are two ways to do this.  The first is modifying the ```<OUTPUT/>``` tag typically found at the end of a script, and the second is by setting flags.

Let's look at a typical usage case for the ```<OUTPUT/>``` tag, first.  Sometimes you may want to use different energy functions during different scoring. For example, you may want to change constraint weights, or to use a lower resolution energy function.  In order to do this, we:

1. Add a named custom scoring function in the ```SCOREFXNS``` section of the XML.

2. Add this to the ```<OUTPUT/>``` tag.

Each custom scorefunction is defined by different sub-tags in the SCOREFXNS section. The format is detailed in the [SCOREFXNS documentation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/RosettaScripts#scorefunctions).

```
<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="molmech" weights="mm_std_fa_elec_dslf_fa13" />
        <ScoreFunction name="r15_cart" weights="ref2015" >
            <Reweight scoretype="pro_close" weight="0.0" />
            <Reweight scoretype="cart_bonded" weight="0.625" />
        </ScoreFunction>
    </SCOREFXNS>
    <RESIDUE_SELECTORS>
    </RESIDUE_SELECTORS>
    <TASKOPERATIONS>
    </TASKOPERATIONS>
    <FILTERS>
    </FILTERS>
    <MOVERS>
    </MOVERS>
    <APPLY_TO_POSE>
    </APPLY_TO_POSE>
    <PROTOCOLS>
    </PROTOCOLS>
    <OUTPUT scorefxn="r15_cart" />
</ROSETTASCRIPTS>
```

The script scorefxn.xml gives and example of defining different scorefunctions. It defines two scorefunctions.  The first one (molmech) is a molecular mechanics scorefunction that is included in the Rosetta database, used as-is, and the second (r15\_cart) is the ref2015 scorefunction modified by changing the weights (coefficients) for certain score terms. (One can also use patch files, or locally-specified weights files; additionally, other scorefunction options can be set, such as soft Lennard-Jones potentials or whatnot.  See the documentation on the ```Set``` tag in the [RosettaScripts documentation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/RosettaScripts) for more on this.)

The molmech scorefunction is never used in this script, which is not a problem -- RosettaScripts does not object to objects that are defined but never used (though the unnecessary allocation of these objects in memory is probably best avoided if one can help it).  The r15\_cart score function *is* used, however, in the OUTPUT tag. This tells RosettaScripts to rescore the output structures with the custom r15\_cart scorefunction, rather than with the default (command line) scorefunction. Run 1ubq.pdb through the script:

```bash
$> cp inputs/scoring.xml .
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol scoring.xml -out:prefix scoring_
```

If you open the scoring_1ubq_0001.pdb output file, you should see that the score table includes columns for the cart\_bonded term, and no pro\_close term.

The above could also be accomplished by passing a custom .wts file to RosettaScripts using the ```-score:weights``` flag at the commandline.

Now let's look at another example of output control at the commandline: we may not want to use PDB output if we're planning to generate very large numbers of structures.  The binary silent file is a proprietary Rosetta format that tis much more compact than a PDB file, and which can store arbitrarily large numbers of structures, avoiding disk space and file count limitations on many file systems.  To produce a silent files for output, let's re-run the command that we just ran, but with an additional flag:

```xml
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol scoring.xml -out:file:silent scoring.silent
```

This time, the output will be a binary silent file.  PDB files can be extracted from binary silent files using the extract\_pdbs application.

## Altering the Pose: Movers

### Minimization

#### Simple Minimization

The core of a RosettaScript XML is the movers. Movers are what will change the structure. Technically, movers are anything that changes the *pose*. While this includes changes to the atomic coordinates, it also includes changes to other features of the pose, including the FoldTree, constraints, sequence, or covalent connectivity. There are certain movers which will change just this auxiliary information, without altering atomic coordinates at all.

* *Minimize the pose before outputting*

As an initial demonstration, we're going to start by writing a script that uses a mover to do gradient-descent energy minimization of the pose. The available movers are listed on the [Rosetta documentation page](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/Movers-RosettaScripts). You can glance through the table of contents for the appropriate section (e.g. "Packing/Minimization") and then look for an appropriate mover (e.g. [MinMover](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/MinMover): Minimizes sidechains and/or backbone).

In each mover documentation page there should be an example tag. Make a copy of scoring.xml called minimization.xml, and then copy that example tag to a place between the ```<MOVERS>``` and ```</MOVERS>``` tags. Again, indentation doesn't matter to Rosetta, although it's easier for *you* to read if things are properly indented.

Depending on the mover, the Mover tags in Rosetta scripts can be configured in one of two ways. They can either take subtags, or they can can have options specified within the tag itself. The MinMover can be configured in both ways. It can take a MoveMap specification as subtags, and other parameters are specified as options (attributes) in the tag itself. What each of the available options means should be described on the mover documentation page.

As mentioned previously, one option that should be in the tag for each mover is the "name" option, with which the user creates a unique handle for referring to that particular instance of the mover (in the PROTOCOLS section, covered below, for example). The value given should be unique to each instance of a mover.  You can have multiple MinMovers as long as their names are different.

The other options in the tag control the mover's behaviour. Most of the options in a tag will have default values associated with them. These are the values which will be used if the option is not provided with the tag. The default values are frequently (though not always) the recommended values for the option, so if you are unsure as to what the option value should be, omitting the option and having it revert to the default value is a good choice. This is what we'll do with the type, tolerance, and max\_iter options in this case. We'll also do this with the MoveMap subtag, leaving it be the default (all atoms move), for now. Other options do not have a default option listed, and if you omit them you will get an error like `Option 'bb' not found in Tag named 'MinMover'`.

Note: Boolean options in the XML can take the same representations of true and false which can be used on the commandline: 1/0, T/F, Y/N, true/false, on/off, etc.

For our example script, we'll make two MinMovers. One we'll call "min\_torsion", which will have the cartesian option set to false (so it will use the default torsional minimization) and will use the molmech scorefunction. The other we'll call "min\_cart", and it will have the cartesian option set to true and use the r15\_cart scorefunction. Both will have bb and chi set to true.

Declaring the movers in the MOVERS section only tells Rosetta that the movers exist and configures their options; however, it doesn't tell Rosetta that they should be applied to the pose (or in what order, or the number of times). The PROTOCOLS section is used to define the sequence of steps that the rosetta\_scripts application will carry out. When RosettaScripts runs on a structure, it will run sequentially through all the entries in the PROTOCOLS section, executing each in order, the output of the previous mover (or filter, as we will see later) becoming the input to the next. In our protocols section we'll add the "min\_cart" mover. Since this is the only mover in the PROTOCOLS section, this is the only mover which will be run. The min\_torsions mover will be defined, but will not be applied to the pose. (The mover can be specified with either the "mover" or "mover\_name" option.)

```
...
    <MOVERS>
        <MinMover name="min_torsion" scorefxn="molmech" chi="true" bb="1" cartesian="F" >
        </MinMover>
        <MinMover name="min_cart" scorefxn="r15_cart" chi="true" bb="1" cartesian="T" >
        </MinMover>
    </MOVERS>
    <APPLY_TO_POSE>
    </APPLY_TO_POSE>
    <PROTOCOLS>
        <Add mover="min_cart" />
    </PROTOCOLS>
...
```

```bash
$> cp inputs/minimize.xml .
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol minimize.xml -out:prefix minimize_
```

Within the tracer output you should see indications that your movers are being used (e.g. "BEGIN MOVER MinMover - min_cart"). Also, if you look at the total scores from the output PDB, you should get much better scores for the minimized 1ubq than the one just rescored with r15\_cart. (about -155 versus +460).

Now let's add the other minimization mover, to demonstrate how movers can be placed in series.  Add the marked line shown below to your script (or use the inputs/minimize2.xml file):

```
...
    <PROTOCOLS>
        <Add mover="min_torsion" /> #Add this line
        <Add mover="min_cart" />
    </PROTOCOLS>
...
```

```bash
$> cp inputs/minimize2.xml .
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol minimize2.xml -out:prefix minimize2_
```

This time, when you run the application, you'll find that the torsion-space minimization is carried out first (using the molecular mechanics scorefunction), and the Cartesian-space minimization is carried out on the output structure from the torsion-space minimization (using the ref2015 scorefunction, modified with the cart\_bonded term turned on and the pro\_close term turned off).  Note that Rosetta does not write out any structures until the end of the protocol.

#### More Advanced Minimization

* *Minimize the pose, controlling the minimization with a MoveMap*.

So far, we have used the MinMover as an example of a generic mover.  The MinMover, however, invokes the Rosetta *minimizer*, which is a fundamental Rosetta algorithm, and so it also serves as a good demonstration of the manner in which we control minimizer behaviour in RosettaScripts.  Movers that use the minimizer typically accept a *MoveMap*.  As discussed in the [minimizer tutorial](../minimization/minimization.md), MoveMaps allow users to set which degrees of freedom (DoFs) are fixed during minimization, and which can be altered by the minimizer.  RosettaScripts provides its own syntax for defining a MoveMap.

> **MoveMaps control the minimizer, and most movers that invoke the minimizer can accept a RosettaScripts-style MoveMap.**

Let's modify the current script to demonstrate how a MoveMap can be set up and used to control the MinMover.  First, let's delete min\_cart from the MOVERS section and its invocation in the PROTOCOLS section, and focus on min\_torsion.  Next, let's add the marked lines, below, to the script (or, alternatively, use the inputs/minimize3.xml file):

```xml
...
    <MOVERS>
        <MinMover name="min_torsion" scorefxn="molmech" chi="true" bb="1" cartesian="F" >
            <MoveMap name="min_torsion_mm">                         # Add this
                <Span begin="1" end="999" chi="false" bb="false" /> # And this
                <Span begin="1" end="50" chi="true" bb="true" />    # And this
                <Span begin="5" end="10" chi="true" bb="false" />   # And this
            </MoveMap>                                              # And this
        </MinMover>
    </MOVERS>
    <APPLY_TO_POSE>
    </APPLY_TO_POSE>
    <PROTOCOLS>
        <Add mover="min_torsion" />
    </PROTOCOLS>
...
```

We're telling the MinMover to make use of a MoveMap that *first* disables sidechain ("chi") and mainchain ("bb") degrees of freedom for all residues, *then* re-enables sidechain and mainchain degrees of freedom for residues 1 through 50, and *then* disables mainchain degrees of freedom for residues 5 through 10.  Note that MoveMaps are not perturbed by poses shorter than the residue ranges in their ```<Span>``` tags, so setting values for residues 1 through 999 is perfectly permissible despite the fact that we're working with a 76-residue structure.  Note also that MoveMaps obey the order of operations given in the tag.  In this example, the final effect is to enable all degrees of freedom for residues 1 through 4, only sidechain degrees of freedom fro residues 5 through 10, all degrees of freedom for residues 11 through 50, and no degrees of freedom for residues 51 through 76.

> **Order of operations matters for MoveMaps.**

We can run this with the following:

```bash
$> cp inputs/minimize3.xml .
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol minimize3.xml -out:prefix minimize3_
```

Practially, it's important to know how to set up MoveMaps because there are many situations in which one may wish to prevent the minimizer from moving parts of a pose.  One example is when designing a binder to a target of known structure: typically, there is little to no advantage to letting the minimizer move the backbone of the target, or sidechains that are far from the binding interface.  Indeed, doing so can result in deceptively low-energy structures with little resemblance to anything physically meaningful.  See the [FastRelax mover's documentation page](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/FastRelaxMover) for full documentation on the MoveMap syntax.

### Packing

#### Repacking sidechains

* *Configure the packer with TaskOperations.*
* *Optimize sidechain conformations using the packer.*

The PackRotamersMover is another commonly-used Rosetta mover.  Because it calls the *packer*, another core Rosetta algorithm (see the [packing tutorial](../Optimizing_Sidechains_The_Packer/Optimizing_Sidechains_The_Packer.md)), the PackRotamersMover is a good mover to use to demonstrate the RosettaScripts interface for controlling the packer.  We do this by defining *TaskOperations*.

> **Just as MoveMaps control the minimizer, TaskOperations control the packer, and movers that invoke the packer will typically accept lists of TaskOperations as inputs.**

Let's create a new skeleton XML, and define the ref2015 scorefunction in it:

```xml
<ROSETTASCRIPTS>
	<SCOREFXNS>
		<ScoreFunction name="r15" weights="ref2015" />
	</SCOREFXNS>
	<RESIDUE_SELECTORS>
	</RESIDUE_SELECTORS>
	<TASKOPERATIONS>
	</TASKOPERATIONS>
	<FILTERS>
	</FILTERS>
	<MOVERS>
	</MOVERS>
	<APPLY_TO_POSE>
	</APPLY_TO_POSE>
	<PROTOCOLS>
	</PROTOCOLS>
	<OUTPUT scorefxn="r15" />
</ROSETTASCRIPTS>
```

In the movers section, let's create a PackRotamersMover.  You can cut-and-paste from the [help page for the mover](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/PackRotamersMover).  Don't forget to add it to the protocols section, as well.  Your MOVERS and PROTOCOLS sections should look something like this:

```xml
...
	<MOVERS>
		<PackRotamersMover name="pack1" scorefxn="r15" task_operations="" />
	</MOVERS>
..
	<PROTOCOLS>
		<Add mover="pack1" />
	</PROTOCOLS>
...
```

Note that, for now, we've left the ```task_operations``` field blank.  Were we to omit this completely (the rosetta\_scripts appliction will throw an error if an option is left blank) and run the script, the PackRotamersMover would call the packer, and the packer would use all rotamers for all 20 canonical amino acids at every position -- that is, it would try to design the entire protein, which is not what we want.

> **The packer's default behaviour is to design with all canonical amino acids at every position.  Preventing design with TaskOperations, or otherwise limiting the behaviour of the packer at some subset of residue positions, is essential for *almost all* usage cases.**

TaskOperations are the means by which the user controls the packer.  They specify which residue to repack and/or design, and how to do it. TaskOperations are defined in the TASKOPERATIONS section of the XML, and as with the movers, the available types are listed on [the corresponding documentation page](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/TaskOperations-RosettaScripts).

In addition to controlling which positions are designed or repacked, TaskOperations also control details about how sidechains are sampled. The default is strictly for on-rotamer sampling, but it's frequently useful to add additional sub-rotameric samples. For example, adding plus or minus one standard deviation around the center of each rotamer bin can help the packer to find better sidechain combinations. The [ExtraRotamersGeneric](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ExtraRotamersGenericOperation) TaskOperation allows you to control the rotamer sampling levels. Generally, adding some additional rotamers to chi1 and chi2 is useful, though the cost is a more complex packing problem and longer convergence time. (There are other ways to control this. For example, the [InitializeFromCommandline](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/InitializeFromCommandlineOperation) task operation allows you to use the -ex1 -ex2 options on the commandline to control rotamer sampling.)

So let's create two TaskOperations.  The first will tell the packer to use only the current amino acid type at each position, and consider only alternative rotamers for that type.  (Technically, this is *disabling* design -- a minor point that will be important later.)  The second will enable some extra rotamers.  In the TASKOPERATIONS section of your script, add a [RestrictToRepacking](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/RestrictToRepackingOperation) TaskOperation and an [ExtraRotamersGeneric](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ExtraRotamersGenericOperation) TaskOperation, giving each a name:

```xml
...
	<TASKOPERATIONS>
		<RestrictToRepacking name="no_design" /> #Note that there are no options except name to set.
		<ExtraRotamersGeneric name="extrachi" ex1="1" ex2="1" ex1_sample_level="1" ex2_sample_level="1" /> #This one allows you to set several options, however.
	</TASKOPERATIONS>
...
```

Down below, in the MOVERS section, let's tell the PackRotamersMover that we created earlier to use these TaskOperations.  Note that we can apply these in any order -- TaskOperations are commutative, which makes them different from MoveMaps.  We'll return to this point later, when we're designing with TaskOperations and ResidueSelectors.

```xml
...
	<MOVERS>
		<PackRotamersMover name="pack1" scorefxn="r15" task_operations="no_design,extrachi" />
	</MOVERS>
...
```

Now let's run this script (or the inputs/repack_only.xml file).  This should generate 2305 rotamers, and take on the order of a second or two to run:

```bash
$> cp inputs/repack_only.xml .
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol repack_only.xml -out:prefix repack_only_
```

If you look at the output, you'll see that the sidechains have been repacked, though in many cases, Rosetta found an optimal rotamer very close to that in the input structure.

#### Advanced packing: Using ResidueSelectors with TaskOperations to redesign the protein core

* *Redesign (*i.e.* find a new sequence for) the ubiquitin core.*
* *Use ResidueSelectors in conjuction with TaskOperations and the PackRotamersMover.*
* *Understand TaskOperation commutativity.*

Let's consider a more complicated (and more realistic) usage case -- one that demonstrates how we can single out subsets of residues in a structure and do different things to different parts of a pose.  Let's find a new sequence for the buried core residues in ubiquitin, while permitting boundary (semi-buried) residues to repack and prohibiting surface residues from moving at all.  We'll also restrict the core to hydrophobic amino acid types.  To do this, we need a way of selecting these layers.  Some of the general TaskOperations are able to select certain residues, but a more flexible choice for selecting certain residues is [ResidueSelectors](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/ResidueSelectors/ResidueSelectors). ResidueSelectors, like their name suggests, are able to specify (select) a particular subset of residues, which can then be used with TaskOperations or other RosettaScripts objects. Unlike TaskOperations, which are strictly one way (you can turn off design, but you can't turn it back on), ResidueSelectors can be combined in various ways to select the particular residue you want.  It's worth taking a moment to comment on the differences between TaskOperations and ResidueSelectors:

| |TaskOperations | ResidueSelectors |
|---|---|---|
| **Intended purpose** | Setting packer behaviours. (*e.g.* Disabling design, limiting allowed residue idenities at certain sequence positions, enabling extra rotamers, telling the packer to include the input rotamer, *etc.*).  Note that, because TaskOperations predate ResidueSelectors, there are some older Rosetta modules that use TaskOperations as a means of selecting residues, though this is being phased out. |  Selecting subsets of residues in a pose based on rules, then passing the subsets as inputs to other Rosetta modules.  |
| **Rule for combining** | Commutativity: applying TaskOperation A, B, and C produces the same effect regardless their order. | Boolean operations: ResidueSelectors produce selections that can be combined to produce the union (OR) or intersection (AND) of the set, or which can be inverted (NOT).  Nested Boolean operations allow very complicated combination rules. |
| **Can be passed to** | Movers that invoke the packer.  (Certain other, older Rosetta modules also accept TaskOperations as a means of selecting residues.  This functionality pre-dates ResidueSelectors, and will at some point be deprecated completely.) | Many movers, filters, and TaskOperations, and even to other ResidueSelectors. |

Let's start by defining three ResidueSelectors to select residues based on burial, in core, boundary, and surface layers.  Of the [available ResidueSelectors](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ResidueSelectors), the [LayerSelector](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ResidueSelectors#residueselectors_conformation-dependent-residue-selectors_layerselector) is the one that will allow us to select residues based on burial (with details of the algorithm available from the [help documentation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ResidueSelectors#residueselectors_conformation-dependent-residue-selectors_layerselector)).  Start a new RosettaScript, define a basic scorefunction, and then define three LayerSelectors in the RESIDUE_SELECTORS section as follows:

```xml
<ROSETTASCRIPTS>
	<SCOREFXNS>
		<ScoreFunction name="r15" weights="ref2015" />
	</SCOREFXNS>
	<RESIDUE_SELECTORS>
		<Layer name="corelayer" select_core="true" select_boundary="false" select_surface="false" core_cutoff="4.0" />
		<Layer name="boundarylayer" select_core="false" select_boundary="true" select_surface="false" core_cutoff="4.0" />
		<Layer name="surfacelayer" select_core="false" select_boundary="false" select_surface="true" core_cutoff="4.0" />
	</RESIDUE_SELECTORS>
	<TASKOPERATIONS>
	</TASKOPERATIONS>
	<FILTERS>
	</FILTERS>
	<MOVERS>
	</MOVERS>
	<APPLY_TO_POSE>
	</APPLY_TO_POSE>
	<PROTOCOLS>
	</PROTOCOLS>
	<OUTPUT scorefxn="r15" />
</ROSETTASCRIPTS>
```

OK, now let's use these to set up some TaskOperations for each layer.  For the core, we want to restrict design to hydrophobic amino acid residue types.  We'll use a [resfile](https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/resfiles) to specify allowed types, and the [ReadResFile TaskOperation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ReadResfileOperation) to apply this restriction to core residues only, using the "corelayer" ResidueSelector to restrict this TaskOperation's action to the core.  Create a new text file (we'll call it core_resfile.txt) and add the following to it:

```
PIKAA FAMILYVW #Pick amino acids PHE, ALA, MET, ILE, LEU, TYR, VAL, or TRP to design with; prohibit all others.
start
```

This is our resfile, indicating that only hydrophobic residues (and alanine) will be allowed.  Note that it only contains a global options line; we're not specifying any per-residue behaviour in this resfile (though that is an option).  Now let's add the ReadResfile TaskOperation.

```xml
...
	<TASKOPERATIONS>
		<ReadResfile name="core_resfile" filename="core_resfile.txt" selector="corelayer" />
	</TASKOPERATIONS>
...
```

Passing "corelayer" with the "selector=" option indicates that, rather than being applied to the whole pose, the effects of the resfile will only be applied to the selected residues.  If this were our only TaskOperation passed to the packer, the overall effect would be to design with all 20 amino acids everywhere *except* in the core, where we would design only with hydrophobic residues.  So now, we need to set the behaviour for boundary and surface layers.  Note, though, that the [RestrictToRepacking TaskOperation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/RestrictToRepackingOperation) that we used earlier takes no options, so there's no direct way to use a ResidueSelector to apply its effect to a subset of residues.  For this reason, we'll use the [OperateOnResidueSubset TaskOperation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/OperateOnResidueSubsetOperation), and two [Residue-Level TaskOperations](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/Residue-Level-TaskOperations): the [RestrictToRepackingRLT](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/Residue-Level-TaskOperations) and the [PreventRepackingRLT](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/Residue-Level-TaskOperations).  As before, we should also define a TaskOperation to allow extra rotamers to be used.

```xml
	<TASKOPERATIONS>
		<ReadResfile name="core_resfile" filename="core_resfile.txt" selector="corelayer" />
		<OperateOnResidueSubset name="restrict_boundary_to_repack" selector="boundarylayer" >
			<RestrictToRepackingRLT />
		</OperateOnResidueSubset>
		<OperateOnResidueSubset name="prevent_surface_from_repackin" selector="surfacelayer" >
			<PreventRepackingRLT />
		</OperateOnResidueSubset>
		<ExtraRotamersGeneric name="extrachi" ex1="1" ex2="1" ex1_sample_level="1" ex2_sample_level="1" />
	</TASKOPERATIONS>
```

The rest is as before: set up a PackRotamersMover, passing the four TaskOperations defined above to it:

```xml
...
	<MOVERS>
		<PackRotamersMover name="pack1" scorefxn="r15" task_operations="core_resfile,prevent_surface_from_repacking,restrict_boundary_to_repack,extrachi" />
	</MOVERS>
...
	<PROTOCOLS>
		<Add mover="pack1" />
	</PROTOCOLS>
...
```

The final file is provided as inputs/design_core.xml.  You can run this with:

```bash
$> cp inputs/design_core.xml .
$> cp inputs/core_resfile.txt .
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol design_core.xml -out:prefix design_core_
```

Open the output.  You might notice that relatively few core residues have changed: Rosetta finds an optimal sequence very similar to the input, as we might expect.  It's also worth noting this line in the output log:

```
core.pack.pack_rotamers: built 4605 rotamers at 29 positions.
```

In comparison, the repacking job that we ran earlier, with no design, generated over 2000 rotamers.  This illustrates an important point: although design is normally a far more computationally expensive task than simple repacking without design, clever use of TaskOperations can reduce the complexity of the problem considerably.  A design job that's set up cleverly can involve comparable computational complexity to a naïve repacking job without design.

> **Pay careful attention to packer setup.  A poorly-conceived packer job can be prohibitively computationally expensive, while a well-designed one can be very quick to execute.**

The task of writing and running a naïve design script, that designs with all 20 amino acid residue types at all positions, to compare rotamers generated and running time, is left as an exercise for the reader.  (If you try this, you'll find that the naïve design run is much, much slower than the one controlled carefully with ResidueSelectors and TaskOperations!)

#### Understanding commutativity of TaskOperations

Let's do an additional thing with the script that we have to illustrate one final point about TaskOperations: let's add one more ReadResfile TaskOperation.  In this second ReadResfile, let's use the PIKAA command to choose a different, but overlapping, set of allowed residue types -- say, PHE, TYR, ASP, GLU, LYS, and ARG.  The resfile (call it core_resfile2.txt) would look like this:

```
start
1 - 76 A PIKAA FYDERK
```

The new ReadResfile TaskOperation, in the TASKOPERATIONS section, would look like this:

```xml
		<ReadResfile name="core_resfile2" filename="core_resfile2.txt" />
```

It should be appended to the list of TaskOperations passed to the PackRotamersMover, like so:

```xml
		<PackRotamersMover name="pack1" scorefxn="r15" task_operations="core_resfile,prevent_surface_from_repacking,restrict_boundary_to_repack,extrachi,core_resfile2" />
```

Run the modified script (or use inputs/design_core2.xml):

```bash
$> cp inputs/design_core2.xml .
$> cp inputs/core_resfile2.txt .
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol design_core2.xml -out:prefix design_core2_
```

This time, if you examine the output, there are several things to note:

1.  The core is now entirely phenylalanine and tryptophan.  This is because canonical residue types can only be turned *off*; once off, they can't be turned back *on*.  This is the AND-commutativity of TaskOperations at work: the packer only designs with a residue type if TaskOperation A *and* TaskOperation B permit it.  Since the core\_resfile TaskOperation prohibits ASP, GLU, LYS, and ARG, and the core\_resfile2 TaskOperation prohibits ALA, MET, ILE, LEU, TYR, and VAL, the only amino acids permitted are TRP and PHE.

2.  Only the core has been designed.  The behaviours of restricting to repacking and preventing repacking override the allowed amino acid types for design, and obey OR-commutativity: if TaskOperation A *or* TaskOperation B indicates that a position should be restricted to repacking or prevented from repacking, then the combination of TaskOperations also results in that residue being restricted to/prevented from repacking.

> **The commutativity of TaskOperations is very important.  Applying A, B, and C is the same as applying C, B, and A.  One must always think carefully about what one is prohibiting or enabling when using combinations of TaskOperations.**

#### Combining ResidueSelectors

As mentioned earlier, ResidueSelectors can be combined with Boolean operations.  This is accomplished with three special ResidueSelectors, called the [AndResidueSelector](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ResidueSelectors#residueselectors_logical-residueselectors_andresidueselector), the [OrResidueSelector](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ResidueSelectors#residueselectors_logical-residueselectors_orresidueselector), and the [NotResidueSelector](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ResidueSelectors#residueselectors_logical-residueselectors_notresidueselector).  The AndResidueSelector and the OrResidueSelector each take as inputs two or more other, previously-defined ResidueSelectors; the selection that they return is the intersection and union of the sets of residues selected by the input ResidueSelectors, respectively (*i.e.* the AndResidueSelector selects a residue if it is selected by input ResidueSelector A *and* input ResidueSelector B, while the OrResidueSelector selects a residue if it is selected by input ResidueSelector A *or* input ResidueSelector B).  The NotResidueSelector inverts a selection, selecting all residues not selected by a single input ResidueSelector.

These Boolean operations allow us to do some very powerful things.  As an example, let's imagine that we were going to modify our first core design script, above, so that now it redesigns the core, but does *not* design or repack existing polar amino acid residues in the core -- we want to preserve those.  To achieve this, we could modify the selector that we pass to the "prevent\_surface\_from\_repacking" TaskOperation, so that it also prevents polar amino acid residues in the core from repacking.

First, we need a ResidueSelector that will select polar amino acid residues.  The ResidueNameSelector will serve nicely for this.  Modify the design_core.xml file (the first script that designed the core, before we did the experiment of adding a second ReadResfile TaskOperation) and add the following to the RESIDUE\_SELECTORS section:

```xml
		<ResidueName name="select_polar" residue_name3="ASP,GLU,LYS,ARG,HIS,SER,THR,ASN,GLN" />
```

Next, let's use an AND selector to select residues that are polar *and* in the core.  We can use the "corelayer" selector that we defined earlier:

```xml
		<And name="polar_and_core" selectors="select_polar,corelayer" />
```

Finally, let's use an OR selector to select residues that are (polar *and* in the core) *or* in the surface layer.  All of these will be restricted to repacking.

```xml
		<Or name="surface_or_buried_polar" selectors="polar_and_core,surfacelayer" />
```

Now, we can pass this selector to the "prevent\_surface\_from\_repacking" TaskOperation, and it will prevent both the surface and the buried polar residues from repacking.  So the definition of the "prevent\_surface\_from\_repacking" TaskOperation changes to:

```xml
		<OperateOnResidueSubset name="prevent_surface_from_repacking" selector="surface_or_buried_polar" >
			<PreventRepackingRLT />
		</OperateOnResidueSubset>
```

The full script is design_core3.xml.  Run it as follows:

```bash
$> cp inputs/core_resfile.txt .
$> cp inputs/design_core3.xml .
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol design_core3.xml -out:prefix design_core3_
```

If you compare the output to the input structure, you'll find that the core has now been redesigned, preserving the buried polar residues' identities and conformations.

#### Summary: Movers

In this sub-section, we have learnt how to set up movers.  In particular, we learnt about movers that call the *minimizer* and which accept *move maps* as inputs, and about movers that call the *packer* and which accept *task operations* as inputs.  Finally, we explored the use of *residue selectors* for defining sets of residues as inputs into other Rosetta modules.

Before we move on to filters, it's worth mentioning that there exist movers that call both the packer and the minimizer over the course of their operation.  These generally accept both move maps and task operations, with the former controlling the minimization steps, and the latter controlling packing steps.  The [FastDesign](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/FastDesignMover) and [FastRelax](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/FastRelaxMover) movers are good examples of this: they both carry out alternating rounds of packing and minimization, gradually ramping the *fa_rep* (Lennard-Jones repulsive) term in the scorefunction from a low value to normal as rounds progress.  It's important to note the division of protocol control between move maps and task operations: if a user wishes to keep a residue completely fixed, he or she must disable movement *both* during packing steps (by preventing repacking with a suitable task operation) *and* during minimization steps (by disabling that residue's degrees of freedom in the move map).

> **Packing and minimization steps in more complicated protocols are controlled by different user interfaces (task operations and move maps, respectively).  Packing and minimization steps must be independently configured for movers and protocols that do both.**

## Filters

* *Filter runs based on a productive conformation (in this case, a salt-bridge)*

Because Rosetta runs are typically stochastic, early stages will often sample conformations which will not be productive. That is, the randomness introduced by initial movers will result in conformations which will never lead to useful final models. To speed up the protocol, it is sometimes helpful to abandon some samples before the final stages of sampling when early stages result in conformations which are known to be unproductive. To facilitate this, RosettaScripts provides Filters, which can stop a job based on measured properties of the protein structure, allowing the rosetta\_scripts application to continue to the next job (*i.e.* the next replicate of the protocol with the current input or the next input structure).

Let's consider the case, now, of repacking just the *surface* (*i.e.* solvent-exposed) residues of ubiqutin, followed by full minimization.  A script to do this might look something like the following:

```xml
<ROSETTASCRIPTS>
	<SCOREFXNS>
		<ScoreFunction name="r15" weights="ref2015" />
	</SCOREFXNS>
	<RESIDUE_SELECTORS>
	</RESIDUE_SELECTORS>
	<TASKOPERATIONS>
		<RestrictToRepacking name="repackonly" />
		<ExtraRotamersGeneric name="extrachi" ex1="1" ex2="1" ex1_sample_level="1" ex2_sample_level="1" />
	</TASKOPERATIONS>
	<FILTERS>
		<AtomicDistance name="salt_bridge" residue1="11A" atomtype1="Nlys" residue2="34A" atomtype2="OOC" distance="3.0" />
	</FILTERS>
	<MOVERS>
		<MinMover name="min" scorefxn="r15" chi="true" bb="true" cartesian="false" />
		<PackRotamersMover name="pack" scorefxn="r15" task_operations="repackonly,extrachi"/>
	</MOVERS>
	<APPLY_TO_POSE>
	</APPLY_TO_POSE>
	<PROTOCOLS>
		<Add mover="pack" />
		<Add filter="salt_bridge" />
		<Add mover="min" />
	</PROTOCOLS>
	<OUTPUT scorefxn="r15" />
</ROSETTASCRIPTS>


```

In this case, we're passing TaskOperations for preventing design and for enabling extra rotamers to a PackRotamersMover.  We also define a MinMover to do minimization.  In the protocols section, we call the PackRotamersMover first, then the MinMover.

In the original structure, there is a salt bridge between K11 and E34, and we probably want to preserve that.  The packer may or may not keep that, though -- sometimes in a packer run, we may not get that. It's the case that if we start with sidechain configurations which are too far apart, minimizing will never pull K11 and E34 back together to re-form the salt bridge.  So if we definitely want the salt bridge in our output structures, the time spent on minimizing the non-salt bridged packing output is effectively wasted.  While this may be seconds in a single run, if we're doing large-scale sampling (say, tens of thousands of trajectories), this could add up to quite a lot of wasted CPU-time.  In many cases, later steps might take minutes or hours, so avoiding unnecessary computation is definitely worthwhile.  Additionally, given that one often manually looks at output structures as a final step, it is good to have a way to reduce the amount of output to a managable number of structures.  In this case, we will use a filter to abandon those jobs that fail to form the salt bridge before we minimize.

> **Filters are important to allow users to abandon non-productive trajectories and to move on to other jobs, to avoid unnecessary computation.**

To enforce the salt bridge in this case, we will filter based off the distance between the two atoms: if they're close enough, we can continue. If they're too far apart, we'll throw out the structure. Skim the [Filters documentation page](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Filters/Filters-RosettaScripts) and look for a filter which might have the appropriate functionality. [AtomicDistance](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Filters/filter_pages/AtomicDistanceFilter) ("Filter based on the distance between two atoms.") looks to be what we want.

As before, copy and paste the example tag from the documentation into the FILTERS section of the XML. As mentioned in the documentation for the filter, you can specify either the specific atom name, or you can specify a Rosetta atom type. If an atom type is specified, then the closest distance for any atom of the relevant type is used. This latter behavior is what we want; we don't care which of the carboxylate oxygens are paired with the lysine side-chain nitrogen. Therefore we can specify the atom types: the "OOC" oxygens from E34 pairing with the "Nlys" nitrogen from K11.

Most filters work by computing some structural metric, and then comparing it to a threshold value to determine if the filter passes or fails. The AtomicDistance filter uses the "distance" options to set the threshold: distances below this pass, distances above fail.

We want to set the distance threshold large enough such that it will pass all the structures which have the salt bridge, but also narrow enough that it will fail the structures which don't have it. (Normally you should err on the side of including too much, as the minimizer may take structures which are slightly outside of the acceptable range and possibly bring them in. However, for this tutorial will use a possibly too narrow distance of 3.0 Ang.)

```xml
...
    <FILTERS>
		<AtomicDistance name="salt_bridge" residue1="11A" atomtype1="Nlys" residue2="34A" atomtype2="OOC" distance="3.0" />
    </FILTERS>
...
```

Again, this only defines the filter. To actually apply it, we have to add it to the protocols section.

```xml
...
    <PROTOCOLS>
        <Add mover="pack" />
        <Add filter="salt_bridge" />
        <Add mover="min" />
    </PROTOCOLS>
...
```

Within the PROTOCOLS section, movers and filters are listed in the order in which they are to be evaluated. That is, the structure will first be packed, then the filter will be applied, and then, if and only if the filter passes, it will be minimized.  If the filter fails, a message is printed to the output log, and the rosetta\_scripts application will continue to the next job.

Let's try this out.  This time, we'll tell the rosetta\_scripts application to repeat the job 100 times with the ```-nstruct 100``` option at the commandline.  We expect that some small fraction of the jobs will succeed and that most will fail to form the salt bridge and will be abandoned.  Note that, by default, Rosetta applications exit with error status if any jobs fail to pass filters.  Since this is not always desireable, there's a way to disable this: the ```-jd2:failed\_job\_exception false``` flag.  (Note that this is essential in MPI mode, since any job failure brings down all processes that are still carrying out jobs.)

```bash
$> cp inputs/filter.xml .
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol filter.xml -out:prefix filter_ -nstruct 100 -jd2:failed_job_exception false
```

Running the above, you'll probably find that about 90% of jobs returned a structure, and 10% failed to pass the filter.  A filter could have any pass rate, though.  This illustrates an important point about filtering: imagine that we had five features that we wanted to filter for, and that each filter passed only 1% of the time.  We would have to do, on average, ten billion samples to obtain one structure.  A better approach is to come up with ways to guide Rosetta to better solutions, increasing the hit rate instead of relying on more sampling and more filtering.  In this case, for example, we could use constraints to guide the packer to form the salt bridge, rather than filtering afterwards.  It's often difficult to come up with simple sampling or scoring biases to use for the desired properties, though, so filtering continues to be a major part of the Rosetta workflow, despite its inefficiency.

> **Filtering abandons non-productive trajectories, but it is more efficient to work to increase the fraction of trajectories that yield productive results than simply to throw away non-productive trajectories.**

Given that job failure is stochastic, this leaves you with fewer output files than you set with -nstruct. You can tell Rosetta to automatically re-run failed jobs by using the `-jd2:ntrials` option. This option sets the number of times each nstruct is retried, if it fails. (It moves on to the next output structure immediately if it was successful.)

### Filters as metric evaluators

In addition to stopping the run, filters can also be used as metric evaluators. For example, we can make filters to compute the heavy atom RMSD of the sidechains for specific residues.  Let's say, for example, that we're interested in the aromatic residues F45 and Y59. From the documentation, it looks like we can use the [SidechainRmsd](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Filters/filter_pages/SidechainRmsdFilter) filter. (We'll use the input structure as the reference pose.) The key to using Filters as metric evaluators instead of as trajectory-stoppers is the "confidence" option for all filters. This tells the filter what random fraction of the time it should act as a filter, and for which it should be just a metric evaluator. The default of "1.0" means always act as a filter. If you set this to "0.0" the filter will never filter, instead it will just act like a metric evaluator, meaning that it reports the value of whatever it calculates, but doesn't ever stop a trajectory based on that value.

Let's add some metric-evaluating filters to the script that we just ran:

```xml
...
	<FILTERS>
		<AtomicDistance name="salt_bridge" residue1="11A" atomtype1="Nlys" residue2="34A" atomtype2="OOC" distance="3.0" />
		<SidechainRmsd name="F45_rmsd" res1_pdb_num="45A" res2_pdb_num="45A" include_backbone="1" confidence="0.0" />
		<SidechainRmsd name="Y59_rmsd" res1_pdb_num="59A" res2_pdb_num="59A" include_backbone="1" confidence="0.0" />
	</FILTERS>
...
	<PROTOCOLS>
		<Add mover="pack" />
		<Add filter="salt_bridge" />
		<Add mover="min" />
		<Add filter="F45_rmsd" />
		<Add filter="Y59_rmsd" />
	</PROTOCOLS>
...
```

Filters used as metric evaluators also need to be added to the PROTOCOLS section. NOTE: While the *filtering* ability of filters take place at their place in PROTOCOLS, the *metric evalution* ability is only applied at the very end of the PROTOCOLS section, to the final, output model.

```bash
$> cp inputs/filter2.xml .
$> rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol filter2.xml -out:prefix filter2_ -nstruct 100 -jd2:failed_job_exception false
```

In addition to printing the results of the metric evaluation to the tracer (output log), the results of the filter will be placed in a column of the scorefile. The name of the column is the same as the name of the filter. Additionally, the values for the filters will be written at the end of the PDB file, after the score table.

## Conclusion

This tutorial was intended to give you a brief introduction to creating an XML protocol. The process we went through is similar to that used by most RosettaScripts developers when writing an XML file from scratch: protocols are built iteratively, starting with a simple protocol and progressively adding different and more complex stages. For each stage, it's important to have an idea about the effect you wish to accomplish, and then to skim the documentation for existing movers/filters/task operations/*etc.* which will accomplish it. This may involve multiple RosettaScripts objects, Rosetta modules that require other Rosetta modules as inputs (*e.g.* movers that require task operations that require residue selectors).

There are, of course, many more RosettaScripts objects than we have discussed, most of which should be covered in the RosettaScripts documentation. There are also additional sections of the XML, which are used for more specialized applications. (For example, ligand docking.)

A final note - even if you can create an XML from scratch, it may be easier not to. If you already have an example XML that does something close to what you want to do, it's probably easier to start with that XML, and alter it to add in the functionality you want.

The hard part is not necessarily in putting together the XML, but in determining the optimal protocol (the logical steps) you should use to accomplish your modeling goals, and then in benchmarking that protocol to make sure it does what you hoped.

## Troubleshooting

RosettaScripts is sensitive to mis-matched tags. If you forget to close a tag, or omit the ending slash on what is supposed to be a standalone tag, RosettaScripts will exit with a (possibly uninformative) error message. If you get something like "Error: Tag::read - parse error", this means there is a syntax error in your XML. The recommended way of debugging it is to make a copy of the script, and progressively comment out or remove portions of the XML file until you get a script that works. (Or at least is able to be parsed.) It is then likely that the source of the error is in the portion of the XML which you commented out or deleted.
