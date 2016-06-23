# Creating protocols with RosettaScripts
======================================

Tutorial by Rocco Moretti (rmorettiase@gmail.com) and Vikram K. Mulligan (vmullig@uw.edu).  Created on 21 June 2016 as part of the 2016 Documentation XRW.

## Goals
--------
At the end of this tutorial, you will understand:
- The RosettaScripts paradigm
- RosettaScripts syntax
- How to control final file output from RosettaScripts
- How to manipulate poses in RosettaScripts using *movers*
	- How to control movers that invoke the minimizer using *MoveMaps*
	- How to control movers that invoke the packer using *TaskOperations*
- How to evaluate pose properties and control protocol flow in RosettaScripts using *filters*
- How to select residues in RosettaScripts using *residue selectors*
- How to nest movers and how to script common loops (*e.g.* Monte Carlo trajectories)
- How to assemble more complicated protocols from simpler building-blocks
- How to control large-scale sampling

## What is RosettaScripts?
-----------------------

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
------------------------

* *Run the simplest possible RosettaScript*

The simplest RosettaScript XML is one which does nothing.  You can obtain a skeleton XML file, which does nothing, in one of two ways.  You go to [the RosettaScripts documentation page](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/RosettaScripts) and find the skeleton XML file there, then copy and paste it into a new file (`nothing.xml`).  You can also generate a skeleton XML file using the rosetta_scripts application:

```bash
$> <path_to_Rosetta_directory>/main/source/bin/rosetta_scripts.default.linuxgccrelease -print_template_script >nothing.xml
```

In the above, the ".default.linuxgccrelease" may need to be changed for your build, operating system, and compiler (*e.g.* ".static.macosclangrelease" for the static build using the clang compiler on the Macintosh operating system).  If you run the above, it will produce the following output:

```
core.init: Rosetta version unknown:979495e360c4960e2a6f41fe9e8bfd5b217e31eb 2016-06-16 21:33:49 -0700 from git@github.com:RosettaCommons/main.git
core.init: command: /home/vikram/rosetta_git/Rosetta/main/source/bin/rosetta_scripts.default.linuxclangrelease -print_template_script
core.init: 'RNG device' seed mode, using '/dev/urandom', seed=34564556 seed_offset=0 real_seed=34564556
core.init.random: RandomGenerator:init: Normal mode, seed=34564556 RG_type=mt19937
core.init: Resolved executable path: /home/vikram/rosetta_git/Rosetta/main/source/build/src/release/linux/3.13/64/x86/clang/3.4-1ubuntu3/default/rosetta_scripts.default.linuxclangrelease
core.init: Looking for database based on location of executable: /home/vikram/rosetta_git/Rosetta/main/database/
apps.public.rosetta_scripts.rosetta_scripts: The -"parser:print_template_script" option was specified.  The app will print a template script and then exit.
apps.public.rosetta_scripts.rosetta_scripts: RosettaScripts script template:

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

This will be written to output.log.  You can delete all lines preceding ```<ROSETTASCRIPTS>``` and following ```</ROSETTASCRIPTS>``` to obtain a minimal template.

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

  Anything not in angle brackets is a comment.

    </MOVERS>
    <APPLY_TO_POSE>
    </APPLY_TO_POSE>
    <PROTOCOLS>
    </PROTOCOLS>
    <OUTPUT />
</ROSETTASCRIPTS>

*(Angle brackets are the greater than/less than signs)
```

The nothing.xml file is also provided in the inputs directory:

```bash
$> cp inputs/nothing.xml .
```

As you haven't further defined any protocol, this XML does nothing to the structure. As a test, let's just run a structure through RosettaScripts with this XML. RosettaScripts takes the standard input and output flags. In addition, the `-parser:protocol` option specifies which XML file to use.

```bash
$> <path_to_Rosetta_directory>/main/source/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol nothing.xml
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
<PackRotamers name=pack1 task_operations=task1,task2,task3 /> #This is allowed
<PackRotamers name=pack2 task_operations=task2, task2, task3 /> #This will be misinterpreted
```

This brings up another RosettaScripts syntax convention: generally, we have blocks that define *types* of objects, and within these blocks, we define individual *instances* of objects of the type, giving each one a unique name.  For example, the ```<MOVERS> ... </MOVERS>``` block is the place to define movers.  Within this, we define specific instances of specific types of movers, and we set options for these movers, including a unique name by which each mover will be addressed at later points in the script.  For example:

```xml
	<MOVERS>  #In this section, movers are defined.
		  #The following is a particular mover of the "PackRotamers" type, which we give the
		  #unique name "pack1".  It takes, as an option, a list of previously-defined
		  #TaskOperation objects (a type of object that will be introduced later in this
		  #tutorial).  We assume that task1, task2, and task3 were defined and given these
		  #unique names prior to this point in the script.
		<PackRotamers name=pack1 task_operations=task1,task2,task3 />
		  #From now on, we can refer to the mover defined above using the unique name "pack1".
	</MOVERS>
```

Looking at the output PDB, the output structure (1ubq\_0001.pdb) should be nearly identical to the input structure. The major difference should be the presence of hydrogens which were not in the input structure. This is *not* something that is specific to RosettaScripts - in general Rosetta will add missing hydrogens and repack sidechain atoms missing in the input PDB.

Additionally, you should see the standard Rosetta score table at the end of the PDB. By default, the structure will be rescored with the default Rosetta score function (talaris2014, as of this writing). This can be controlled by the ```-score:weights``` command line option. 

## Controlling RosettaScripts File Output
--------------

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
        <t13 weights="talaris2013" />
        <t14_cart weights="talaris2014" >
            <Reweight scoretype="pro_close" weight="0.0" />
            <Reweight scoretype="cart_bonded" weight="0.625" />
        </t14_cart>
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
    <OUTPUT scorefxn="t14_cart" />
</ROSETTASCRIPTS>
``` 

The script scorefxn.xml gives and example of defining different scorefunctions. It defines two scorefunctions.  The first one (t13) is simply the talaris2013 weights used as-is, and the second is the talaris2014 weights modified by changing the weights (coefficients) for certain score terms. (One can also use patch files, or locally-specified weights file; additionally, other scorefunction options can be set, such as soft Lennard-Jones potentials or whatnot.  See the documentation on the ```Set``` tag in the [RosettaScripts documentation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/RosettaScripts) for more on this.)

The t13 scorefunction is never used in this script, which is not a problem -- RosettaScripts does not object to objects that are defined but never used (though the unnecessary allocation of these objects in memory is probably best avoided if one can help it).  The t14\_cart score function *is* used, however, in the OUTPUT tag. This tells RosettaScripts to rescore the output structures with the custom t14\_cart score function, rather than with the default (command line) scorefunction. Run 1ubq.pdb through the script: 

```bash
$> cp inputs/scoring.xml .
$> <path_to_Rosetta_directory>/main/source/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol scoring.xml -out:prefix scoring_
```

If you open the scoring_1ubq_0001.pdb output file, you should see that the score table includes columns for the cart\_bonded term, and no pro\_close term.

The above could also be accomplished by passing a custom .wts file to RosettaScripts using the ```-score:weights``` flag at the commandline.

Now let's look at another example of output control at the commandline: we may not want to use PDB output if we're planning to generate very large numbers of structures.  The binary silent file is a proprietary Rosetta format that tis much more compact than a PDB file, and which can store arbitrarily large numbers of structures, avoiding disk space and file count limitations on many file systems.  To produce a silent files for output, let's re-run the command that we just ran, but with an additional flag:

```xml
$> <path_to_Rosetta_directory>/main/source/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol scoring.xml -out:file:silent scoring.silent
```

This time, the output will be a binary silent file.  PDB files can be extracted from binary silent files using the extract\_pdbs application.

## Altering the Pose: Movers
----------------------------

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

For our example script, we'll make two MinMovers. One we'll call "min\_torsion", which will have the cartesian option set to false (so it will use the default torsional minimization) and will use the t13 scorefunction. The other we'll call "min\_cart", and it will have the cartesian option set to true and use the t14\_cart scorefunction. Both will have bb and chi set to true.

Declaring the movers in the MOVERS section only tells Rosetta that the movers exist and configures their options; however, it doesn't tell Rosetta that they should be applied to the pose (or in what order, or the number of times). The PROTOCOLS section is used to define the sequence of steps that the rosetta\_scripts application will carry out. When RosettaScripts runs on a structure, it will run sequentially through all the entries in the PROTOCOLS section, executing each in order, the output of the previous mover (or filter, as we will see later) becoming the input to the next. In our protocols section we'll add the "min\_cart" mover. Since this is the only mover in the PROTOCOLS section, this is the only mover which will be run. The min\_torsions mover will be defined, but will not be applied to the pose. (The mover can be specified with either the "mover" or "mover\_name" option.) 

```
...
    <MOVERS>
        <MinMover name="min_torsion" scorefxn="t13" chi="true" bb="1" cartesian="F" >
        </MinMover>
        <MinMover name="min_cart" scorefxn="t14_cart" chi="true" bb="1" cartesian="T" >
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
$> <path_to_Rosetta_directory>/main/source/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol minimize.xml -out:prefix minimize_
```

Within the tracer output you should see indications that your movers are being used (e.g. "BEGIN MOVER MinMover - min_cart"). Also, if you look at the total scores from the output PDB, you should get much better scores for the minimized 1ubq than the one just rescored with t14_cart. (about -155 versus +460).

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
$> <path_to_Rosetta_directory>/main/source/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol minimize2.xml -out:prefix minimize2_
```

This time, when you run the application, you'll find that the torsion-space minimization is carried out first (using the talaris\_2013 scorefunction), and the Cartesian-space minimization is carried out on the output structure from the torsion-space minimization (using the talaris\_2014 scorefunction, modified with the cart\_bonded term turned on and the pro\_close term turned off).  Note that Rosetta does not write out any structures until the end of the protocol.

#### More Advanced Minimization

* *Minimize the pose, controlling the minimization with a MoveMap*.

So far, we have used the MinMover as an example of a generic mover.  The MinMover, however, invokes the Rosetta *minimizer*, which is a fundamental Rosetta algorithm, and so it also serves as a good demonstration of the manner in which we control minimizer behaviour in RosettaScripts.  Movers that use the minimizer typically accept a *MoveMap*.  As discussed in the [minimizer tutorial](../minimization/minimization.md), MoveMaps allow users to set which degrees of freedom (DoFs) are fixed during minimization, and which can be altered by the minimizer.  RosettaScripts provides its own syntax for defining a MoveMap.

> **MoveMaps control the minimizer, and most movers that invoke the minimizer can accept a RosettaScripts-style MoveMap.**

Let's modify the current script to demonstrate how a MoveMap can be set up and used to control the MinMover.  First, let's delete min\_cart from the MOVERS section and its invocation in the PROTOCOLS section, and focus on min\_torsion.  Next, let's add the marked lines, below, to the script (or, alternatively, use the inputs/minimize3.xml file):

```xml
...
    <MOVERS>
        <MinMover name="min_torsion" scorefxn="t13" chi="true" bb="1" cartesian="F" >
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

We're telling the MinMover to make use of a MoveMap that *first* disables side-chain ("chi") and main-chain ("bb") degrees of freedom for all residues, *then* re-enables side-chain and main-chain degrees of freedom for residues 1 through 50, and *then* disables main-chain degrees of freedom for residues 5 through 10.  Note that MoveMaps are not perturbed by poses shorter than the residue ranges in their ```<Span>``` tags, so setting values for residues 1 through 999 is perfectly permissible despite the fact that we're working with a 76-residue structure.  Note also that MoveMaps obey the order of operations given in the tag.  In this example, the final effect is to enable all degrees of freedom for residues 1 through 4, only side-chain degrees of freedom fro residues 5 through 10, all degrees of freedom for residues 11 through 50, and no degrees of freedom for residues 51 through 76.

> **Order of operations matters for MoveMaps.**

We can run this with the following:

```bash
$> cp inputs/minimize3.xml .
$> <path_to_Rosetta_directory>/main/source/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol minimize3.xml -out:prefix minimize3_
```

Practially, it's important to know how to set up MoveMaps because there are many situations in which one may wish to prevent the minimizer from moving parts of a pose.  One example is when designing a binder to a target of known structure: typically, there is little to no advantage to letting the minimizer move the backbone of the target, or side-chains that are far from the binding interface.  Indeed, doing so can result in deceptively low-energy structures with little resemblance to anything physically meaningful.  See the [FastRelax mover's documentation page](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/FastRelaxMover) for full documentation on the MoveMap syntax.

### Packing

#### Repacking side-chains

* *Configure the packer with TaskOperations.*
* *Optimize side-chain conformations using the packer.*

The PackRotamersMover is another commonly-used Rosetta mover.  Because it calls the *packer*, another core Rosetta algorithm (see the [packing tutorial](../Optimizing_Sidechains_The_Packer/Optimizing_Sidechains_The_Packer.md)), the PackRotamersMover is a good mover to use to demonstrate the RosettaScripts interface for controlling the packer.  We do this by defining *TaskOperations*.

> **Just as MoveMaps control the minimizer, TaskOperations control the packer, and movers that invoke the packer will typically accept lists of TaskOperations as inputs.**

Let's create a new skeleton XML, and define the talaris2014 scorefunction in it:

```xml
<ROSETTASCRIPTS>
	<SCOREFXNS>
		<t14 weights="talaris2014" />
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
	<OUTPUT scorefxn="t14" />
</ROSETTASCRIPTS>
```

In the movers section, let's create a PackRotamersMover.  You can cut-and-paste from the [help page for the mover](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/PackRotamersMover).  Don't forget to add it to the protocols section, as well.  Your MOVERS and PROTOCOLS sections should look something like this:

```xml
...
	<MOVERS>
		<PackRotamersMover name="pack1" scorefxn=t14 task_operations= />
	</MOVERS>
..
	<PROTOCOLS>
		<Add mover="pack1" />
	</PROTOCOLS>
...
```

Note that, for now, we've left the ```task_operations``` field blank.  Were we to omit this and run the script, the PackRotamersMover would call the packer, and the packer would use all rotamers for all 20 canonical amino acids at every position -- that is, it would try to design the entire protein, which is not what we want.

> **The packer's default behaviour is to design with all canonical amino acids at every position.  Preventing design with TaskOperations, or otherwise limiting the behaviour of the packer at some subset of residue positions, is essential for *almost all* usage cases.**

TaskOperations are the means by which the user controls the packer.  They specify which residue to repack and/or design, and how to do it. TaskOperations are defined in the TASKOPERATIONS section of the XML, and as with the movers, the available types are listed on [the corresponding documentation page](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/TaskOperations-RosettaScripts).

In addition to controlling which positions are designed or repacked, TaskOperations also control details about how sidechains are sampled. The default is strictly for on-rotamer sampling, but it's frequently useful to add additional sub-rotameric samples. For example, adding plus or minus one standard deviation around the center of each rotamer bin can help the packer to find better side-chain combinations. The [ExtraRotamersGeneric](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ExtraRotamersGenericOperation) TaskOperation allows you to control the rotamer sampling levels. Generally, adding some additional rotamers to chi1 and chi2 is useful, though the cost is a more complex packing problem and longer convergence time. (There are other ways to control this. For example, the [InitializeFromCommandline](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/InitializeFromCommandlineOperation) task operation allows you to use the -ex1 -ex2 options on the commandline to control rotamer sampling.)

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
		<PackRotamersMover name="pack1" scorefxn="t14" task_operations="no_design,extrachi" />
	</MOVERS>
...
```

Now let's run this script (or the inputs/repack_only.xml file).  This should generate 6134 rotamers, and take on the order of 10 seconds to run:

```bash
$> cp inputs/repack_only.xml .
$> <path_to_Rosetta_directory>/main/source/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol repack_only.xml -out:prefix repack_only_
```

If you look at the output, you'll see that the side-chains have been repacked, though in many cases, Rosetta found an optimal rotamer very close to that in the input structure.

#### Advanced packing: Using ResidueSelectors with TaskOperations to redesign the protein core

* *Redesign (*i.e.* find a new sequence for) the ubiquitin core.*
* *Use ResidueSelectors in conjuction with TaskOperations and the PackRotamersMover.*
* *Understand TaskOperation commutativity.*

Let's consider a more complicated (and more realistic) usage case -- one that demonstrates how we can single out subsets of residues in a structure and do different things to different parts of a pose.  Let's find a new sequence for the buried core residues in ubiquitin, while permitting boundary (semi-buried) residues to repack and prohibiting surface residues from moving at all.  We'll also restrict the core to hydrophobic amino acid types.  To do this, we need a way of selecting these layers.  Some of the general TaskOperations are able to select certain residues, but a more flexible choice for selecting certain residues is [ResidueSelectors](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ResidueSelectors). ResidueSelectors, like their name suggests, are able to specify (select) a particular subset of residues, which can then be used with TaskOperations or other RosettaScripts objects. Unlike TaskOperations, which are strictly one way (you can turn off design, but you can't turn it back on), ResidueSelectors can be combined in various ways to select the particular residue you want.  It's worth taking a moment to comment on the differences between TaskOperations and ResidueSelectors:

| |TaskOperations | ResidueSelectors |
|---|---|---|
| **Intended purpose** | Setting packer behaviours. (*e.g.* Disabling design, limiting allowed residue idenities at certain sequence positions, enabling extra rotamers, telling the packer to include the input rotamer, *etc.*).  Note that, because TaskOperations predate ResidueSelectors, there are some older Rosetta modules that use TaskOperations as a means of selecting residues, though this is being phased out. |  Selecting subsets of residues in a pose based on rules, then passing the subsets as inputs to other Rosetta modules.  |
| **Rule for combining** | Commutativity: applying TaskOperation A, B, and C produces the same effect regardless their order. | Boolean operations: ResidueSelectors produce selections that can be combined to produce the union (OR) or intersection (AND) of the set, or which can be inverted (NOT).  Nested Boolean operations allow very complicated combination rules. |
| **Can be passed to** | Movers that invoke the packer | Many movers, filters, and TaskOperations, and even to other ResidueSelectors. |

Let's start by defining three ResidueSelectors to select residues based on burial, in core, boundary, and surface layers.  Of the [available ResidueSelectors](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ResidueSelectors), the [LayerSelector](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ResidueSelectors#residueselectors_conformation-dependent-residue-selectors_layerselector) is the one that will allow us to select residues based on burial (with details of the algorithm available from the [help documentation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ResidueSelectors#residueselectors_conformation-dependent-residue-selectors_layerselector)).  Start a new RosettaScript, define a basic scorefunction, and then define three LayerSelectors in the RESIDUE_SELECTORS section as follows:

```xml
<ROSETTASCRIPTS>
	<SCOREFXNS>
		<t14 weights="talaris2014" />
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
	<OUTPUT scorefxn="t14" />
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
		</OpearateOnResidueSubset>
		<OperateOnResidueSubset name="prevent_surface_from_repackin" selector="surfacelayer" >
			<PreventRepackingRLT />
		</OpearateOnResidueSubset>
		<ExtraRotamersGeneric name="extrachi" ex1="1" ex2="1" ex1_sample_level="1" ex2_sample_level="1" />
	</TASKOPERATIONS>
```

The rest is as before: set up a PackRotamersMover, passing the four TaskOperations defined above to it:

```xml
...
	<MOVERS>
		<PackRotamersMover name="pack1" scorefxn="t14" task_operations="core_resfile,prevent_surface_from_repacking,restrict_boundary_to_repack,extrachi" />
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
$> cp core_resfile.txt .
$> <path_to_Rosetta_directory>/main/source/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol design_core.xml -out:prefix design_core_
```

Open the output.  You might notice that relatively few core residues have changed: Rosetta finds an optimal sequence very similar to the input, as we might expect.  It's also worth noting this line in the output log:

```
core.pack.pack_rotamers: built 8928 rotamers at 50 positions.
```

In comparison, the repacking job that we ran earlier, with no design, generated over 6000 rotamers.  This illustrates an important point: although design is normally a far more computationally expensive task than simple repacking without design, clever use of TaskOperations can reduce the complexity of the problem considerably.  A design job that's set up cleverly can involve comparable computational complexity to a naïve repacking job without design.

> **Pay careful attention to packer setup.  A poorly-conceived packer job can be prohibitively computationally expensive, while a well-designed one can be very quick to execute.**

#### Understanding commutativity of TaskOperations

Let's do one more thing to illustrate one final point about TaskOperations: let's add one more ReadResfile TaskOperation.  In this second ReadResfile, let's use the PIKAA command to choose a different, but overlapping, set of allowed residue types -- say, PHE, TYR, ASP, GLU, LYS, and ARG.  The resfile (call it core_resfile2.txt) would look like this:

```
start
1 - 76 A PIKAA FWDERK
```

The new ReadResfile TaskOperation, in the TASKOPERATIONS section, would look like this:

```xml
		<ReadResfile name="core_resfile2" filename="core_resfile2.txt" />
```

It should be appended to the list of TaskOperations passed to the PackRotamersMover, like so:

```xml
		<PackRotamersMover name="pack1" scorefxn="t14" task_operations="core_resfile,prevent_surface_from_repacking,restrict_boundary_to_repack,extrachi,core_resfile2" />
```

Run the modified script (or use inputs/design_core2.xml):

```bash
$> cp inputs/design_core2.xml .
$> cp core_resfile2.txt .
$> <path_to_Rosetta_directory>/main/source/bin/rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol design_core2.xml -out:prefix design_core2_
```

This time, if you examine the output, there are several things to note:
1.  The core is now entirely phenylalanine and tryptophan.  This is because canonical residue types can only be turned *off*; once off, they can't be turned back *on*.  This is the AND-commutativity of TaskOperations at work: the packer only designs with a residue type if TaskOperation A *and* TaskOperation B permit it.  Since the core\_resfile TaskOperation prohibits ASP, GLU, LYS, and ARG, and the core\_resfile2 TaskOperation prohibits ALA, MET, ILE, LEU, TYR, and VAL, the only amino acids permitted are TRP and PHE.
2.  Only the core has been designed.  The behaviours of restricting to repacking and preventing repacking override the allowed amino acid types for design, and obey OR-commutativity: if TaskOperation A *or* TaskOperation B indicates that a position should be restricted to repacking or prevented from repacking, then the combination of TaskOperations also results in that residue being restricted to/prevented from repacking.

> **The commutativity of TaskOperations is very important.  Applying A, B, and C is the same as applying C, B, and A.  One must always think carefully about what one is prohibiting or enabling when using combinations of TaskOperations.**

CONTINUE HERE

# Filters
---------

* *Filter runs based on a productive conformation (e.g. a salt-bridge)*

Because Rosetta runs are typically stochastic, early stages will often sample conformations which will not be productive. That is, the randomness introduced by initial movers will result in conformations which will never lead to useful final models. To speed up the protocol, it is sometimes helpful to abandon some samples before the final stages of sampling when early stages result in conformations which are known to be unproductive. To facilitate this, RosettaScripts provides Filters, which can stop a protocol based on measured properties of the protein structure.

In our sample packing run, we sometimes get a salt bridge between R54 and D58, but frequently we don't. It's the case that if we start with sidechain configurations which are too far apart, the minimizer will never build the salt bridge. So if we definitely want the salt bridge in our output structures, the time spent on minimizing the non-salt bridged packing output is effectively wasted.

(Note that in a real run we might be better off using constraints to bias the score function used in packing such that the desired hydrogen bond receives a bonus, rather than filtering afterwards. In general, it's normally more efficient to bias sampling and scoring during structure generation, rather than attempt to filter out structures with bad geometries later. However, it's often difficult to come up with simple sampling or scoring biases to use for the desired properties, so filtering is the most straightforward way to do it.)

To enforce the salt bridge, we want to filter based off the distance between the two atoms: if they're close enough, we can continue. If they're too far apart, we'll throw out the structure. Scan the [Filters documentation page](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Filters/Filters-RosettaScripts) and look for a filter which might have the appropriate functionality. [AtomicDistance](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Filters/filter_pages/AtomicDistanceFilter) ("Filter based on the distance between two atoms.") looks to be what we want.

As before, copy and paste the example tag from the documentation into the FILTERS section of the XML. As mentioned in the documentation for the filter, you can specify either the specific atom name, or you can specify a Rosetta atom type. If an atom type is specified, then the closest distance for any atom of the relevant type is used. This latter behavior is what we want; we don't care which of the carboxylate oxygens are paired with which of the guanidinium nitrogens. Therefore we can specify the atom types: the "OOC" oxygens from D58 pairing with the "Narg" nitrogens from R54.

Most filters work by computing some structural metric, and then comparing it to a threshold value to determine if the filter passes or fails. The AtomicDistance filter uses the "distance" options to set the threshold: distances below this pass, distances above fail. 
We want to set the distance threshold large enough such that it will pass all the structures which have the salt bridge, but also narrow enough that it will fail the structures which don't have it. (Normally you should err on the side of including too much, as the minimizer may take structures which are slightly outside of the acceptable range and possibly bring them in. However, for this tutorial will use a possibly too narrow distance of 3.0 Ang.)

```
    <FILTERS>
        <AtomicDistance name="salt_bridge" residue1="54A" atomtype1="Narg" residue2="58A" atomtype2="OOC" distance="3.0" />
    </FILTERS>
```

Again, this only defines the filter. To actually apply it, we have to add it to the protocols section.

```
    <PROTOCOLS>
        <Add mover="pack" />
        <Add filter="salt_bridge" />
        <Add mover="min_cart" />
    </PROTOCOLS>
```

Within the PROTOCOLS section, things are provided in the order they are evaluated. That is, the structure will first be packed, then the filter will be applied, and then it will be minimized.

### Filters as metric evaluators

In addition to stopping the run, filters can also be used as metric evaluators. For example, we can make filters to compute the heavy atom RMSD of the sidechains for specific residues (e.g. F45 and Y59). From the documentation, it looks like we can use the [SidechainRmsd](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Filters/filter_pages/SidechainRmsdFilter) filter. (We'll use the input structure as the reference pose.) The key to using Filters as metric evaluators instead of filters is the "confidence" option for all filters. This tells the filter what random fraction of the time it should act as a filter, and for which it should be just a metric evaluator. The default of "1.0" means always act as a filter. If you set this to "0.0" the filter will never filter, instead it will just act like a metric evaluator.

```
    <FILTERS>
        <SidechainRmsd name="F45_rmsd" res1_pdb_num="45A" res2_pdb_num="45A" include_backbone="1" confidence="0.0" />
        <SidechainRmsd name="Y59_rmsd" res1_pdb_num="59A" res2_pdb_num="59A" include_backbone="1" confidence="0.0" />
    </FILTERS>
```

Filters used as metric evaluators also need to be added to the PROTOCOLS section. NOTE: While the *filtering* ability of filters take place at their place in PROTOCOLS, the *metric evalution* ability is only applied at the very end of the PROTOCOLS section, to the final, output model.

	$> cp inputs/filter.xml .
	$> rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol filter.xml -out:prefix filter_ -nstruct 5

This should result in a fair number of "failed" jobs. (e.g. "ERROR: Exception caught by rosetta_scripts application:3 jobs failed; check output for error messages") This is because the filter will recognize that a large number of generated structures don't match the desired parameters, and will cancel the job.

Given that job failure is stochastic, this will leave you with fewer output files than you set with -nstruct. You can tell Rosetta to automatically re-run failed jobs by using the `-jd2:ntrials` option. This option sets the number of times each nstruct is retried, if it fails. (It moves on to the next output structure immediately if it was successful.)

	$> rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol filter.xml -out:prefix filter2_ -nstruct 5 -jd2:ntrials 10

This should give you five output structures, even if some tries through failed (for example, you'll get messages in the tracer like "filter2_1ubq_0001 reported failure and will retry" and "5 jobs considered, 7 jobs attempted".

In addition to printing the results of the metric evaluation to the tracer, the results of the filter will be placed in a column of the scorefile. The name of the column is the same as the name of the filter. Additionally, the values for the filters will be output to the end of the PDB, after the score table.

Nesting movers
--------------

* *Loop over sidechain optimization until the score doesn't improve.*

One of the more powerful parts of RosettaScripts is the ability to combine individual components in flexible ways. You saw some of this above, where we used ResidueSelectors and TaskOperations as parameters to the PackRotamers mover. There are also certain movers which can take other movers as parameters. This can be used to implement looping.

For our example protocol, we'll add the [RotamerTrialsMinMover](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/RotamerTrialsMinMover), which loops through each residue position, exhaustively testing each position to see if a rotamer substitution will improve things. However, as the ideal sidechain conformation depends on the other sidechains, so the results of a RotamerTrialsMinMover depends on the (random) order in which the sidechains are tested. To make sure we get the best score we possibly can, we're going to repeat the RotamerTrialsMinMover until the score function doesn't improve. 

To do this, we'll use the [IteratedConvergence](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/IteratedConvergenceMover) mover. This mover is a "meta mover" in the sense that it doesn't change the pose itself, but takes as a parameter a Mover which does. It also takes a filter, which is used as a metric evaluator. The IteratedConvergence mover repeatedly applies the given mover, and after each application will call the metric evaluation property of the filter. If the mover keeps improving the score, the IteratedConvergence mover will keep calling the mover. If not, it will stop and return the updated pose.

For the filter, we'll use the [ScoreType](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Filters/filter_pages/ScoreTypeFilter) filter to get the total score of the pose. Since the IteratedConvergence mover only uses this as a metric evaluator, we don't need to worry too much about the threshold or the confidence setting.
 
```
    <FILTERS>
        <ScoreType name="total_score" scorefxn="t14_cart" score_type="total_score" threshold="0"/>
    </FILTERS>
    <MOVERS>
        <RotamerTrialsMinMover name="rtmm" scorefxn="t14_cart" task_operations="repackonly,extrachi,nopack_F45_Y59" />
        <IteratedConvergence name="rotopt" mover="rtmm" filter="total_score" delta="0.1" cycles="1" />
    </MOVERS>
```

Note that when you nest movers/filters/etc. the definition of the sub-mover/filter/etc. must come before the point of use. (Otherwise the order of definition shouldn't matter.) This might involve you making multiple MOVERS/FILTERS/etc. section.

	$> cp inputs/pack_opt.xml .
	$> rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol pack_opt.xml -out:prefix packopt_ -nstruct 2 -jd2:ntrials 10

Looking at the tracer output, you should be able to see the application of the IteratedConvergence, and how the RotamerTrialsMinMover is repeated multiple times.

Variable substition: adding variables to scripts
------------------------------------------------

Sometimes in a RosettaScripts protocol, you want to vary the options given to the tags. For example, if you wish to do a series of runs, with changes at different residues. The naive way of doing this is to make separate XMLs, one for each variant of the option. If you have a large number of variants, this may be less than ideal.

To accomodate this sort of protocol, RosettaScripts has [variable substition](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/RosettaScripts#options-available-in-the-xml-protocol-file_variable-substitution). Within the script you add "%%var_name%%" instead of the option value, and then use the "-parser:script_vars" command line option to set it from the command line.

(NOTE: The variable substitution is only intended for substituting individual options in a tag. Don't try to use it to substitute entire sub-tags.)

For our sample protocol, let's run a mutational scan. There are several movers which can do mutational scanning, but for the purposes of introducing the script_vars functionality, let's use [MutateResidue](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/MutateResidueMover). Also, to keep the runtime short, let's disable the rotamer optimization.

```
    <MOVERS>
        <MutateResidue name="mutate" target="%%position%%" new_res="%%res%% />  
    </MOVERS>
```

To run, we need to then pass something like "-parser:script_vars position=14A new_res=ALA" on the commandline.

	$> cp inputs/mut_scan.xml .
	$> rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol mut_scan.xml -out:prefix V5W_ -nstruct 1 -parser:script_vars position=5A res=TRP -jd2:ntrials 10
	$> rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol mut_scan.xml -out:prefix L43W_ -nstruct 1 -parser:script_vars position=43A res=TRP -jd2:ntrials 10
	$> rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol mut_scan.xml -out:prefix L56W_ -nstruct 1 -parser:script_vars position=56A res=TRP -jd2:ntrials 10
	$> rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol mut_scan.xml -out:prefix L67W_ -nstruct 1 -parser:script_vars position=67A res=TRP -jd2:ntrials 10

These commands should produce a tryptophan scan of a selection of residues in the core of the protein. (Open up the structures in PyMol or the equivalent and compare.

If you wish to do a more thorough scan, either of more positions or of more residue identities, you can easily automate running of the scan by using shell scripting.

## ResidueSelectors and TaskOperations
-----------------------------------


### ResidueSelectors

Looking at available ResidueSelectors, the [ResidueIndexSelector](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ResidueSelectors#residueselectors_conformation-independent-residue-selectors_residueindexselector) looks to be the one we want, when we want to select particular residues. We can specify which residues to use with a comma separated list. Note that we can either use *Pose numbering* (numbers without a chain letter) or *PDB numbering* (with a chain letter). If the PDB hasn't been renumbered to match Pose numbering, these will be different. 

```
    <RESIDUE_SELECTORS>
        <Index name="key_residues" resnums="45A,59A"/>
    </RESIDUE_SELECTORS>
```

We're also only interested in keeping these residues from repacking if they stay as phenylalanine or tyrosine. Let's add a ResidueSelector which selects only phenylalanine or tyrosine residues. Scanning through the documentation, the [ResidueNameSelector](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ResidueSelectors#residueselectors_conformation-independent-residue-selectors_residuenameselector) is a good candidate. We can specify selecting residues with the appropriate three letter codes. 

Note that ResidueSelectors select based on the properties of the structure at the time which they are applied, not the input pose. This means that positions which start as phenylalanine but mutate to a different amino acid will be selected by a PHE ResidueNameSelector used before mutation, but won't be selected if the ResidueSelector is used after the mutation. 

```
    <RESIDUE_SELECTORS>
        <Index name="key_residues" resnums="45A,59A"/>
	<ResidueName name="phe_tyr" residue_name3="PHE,TYR" />
    </RESIDUE_SELECTORS>
```

We want to combine these two selectors. We want selectors which are both PHE or TYR *and* are at position 45A or 59A. To combine ResidueSelectors, we can use one of the "logical" ResidueSelectors. Specifically, we want the [AndResidueSelector](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/ResidueSelectors#residueselectors_logical-residueselectors_andresidueselector), which will select only those residues which are in all of the sub-residue selectors. (Both key_residues *and* phe_tyr.) Residues selected by one *or* the other (but not both) won't be selected by the combined selector. (For example, F4 is selected by phe_tyr, but not by key_residues, so it won't be selected by the joint ResidueSelector.

There's two ways to specify which residue selectors to combine: we can either give a previously defined ResidueSelector by name in the tag, or define new ones as subtags. If we go the subtag route, we don't need to give the ResidueSelectors names. (This is why the "name" option is not listed in the ResidueSelector tag example in the documentation page.)

```
    <RESIDUE_SELECTORS>
	<And name="F45_Y59" >
            <Index resnums="45A,59A" />
	    <ResidueName residue_name3="PHE,TYR" />
        </And>
    </RESIDUE_SELECTORS>
```   

### ResidueSelectors and TaskOperations

The previous section only defined the residue selector - it didn't specify how it was to be used. In our case, we want to use the ResidueSelector to turn off packing to the given residues. Controlling packing is done with TaskOperations, so we need a TaskOperation which can use ResidueSelectors. Looking at the [TaskOperation documentation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/TaskOperations-RosettaScripts), it looks like the [OperateOnResidueSubset](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/OperateOnResidueSubsetOperation) TaskOperation is what we want. This takes a ResidueSelector to define which residues it operates over, and a [ResidueLevelTaskOperation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/TaskOperations/taskoperations_pages/Residue-Level-TaskOperations) to specify what to do with those residues. In our case, we want to use the previously named ResidueSelector, and prevent repacking (so the PreventRepackingRLT ResidueLevelTask).

```
    <TASKOPERATIONS>
        <OperateOnResidueSubset name="nopack_F45_Y59" selector="F45_Y59" >
	    <PreventRepackingRLT/>
        </OperateOnResidueSubset>
    </TASKOPERATIONS>
```  

### TaskOperations and Movers

Again, that just defined the TaskOperation. To use it, we need to pass it to another object (e.g. a mover) which will apply it. If you remember, we were defining a PackRotamersMover. This takes a comma-separated list of TaskOperation names. These TaskOperations will be combined in the standard restrictive fashion. That is, you'll start with all positions set to design to all canonical amino acids, you turn off design and repacking at particular positions, and once design or repacking is turned off, it stays turned off and can't get turned back on.

So we're going to combine the two task operations we defined earlier. The repackonly TaskOperation will turn off design at all positions (including F45 and Y59) and the nopack_F45_Y59 operation, which will turn off repacking (and design) specifically at the selected residues.

```
    <MOVERS>
        <PackRotamersMover name="pack" scorefxn="t13" task_operations="repackonly,extrachi,nopack_F45_Y59"/>
    </MOVERS>
```

As before, putting the tag in the MOVERS section only defines the mover - in order to actually apply the mover, we need to put it in the PROTOCOLS section. The order in which we place the movers in the PROTOCOLS section matters, as the output of one mover will be used as the input to the next. So there's a difference between packing and then minimizing and minimizing and then packing. Typically, you would want to pack and then minimize, as packing does a courser, wider sampling, while minimization is a more local refinement.

```
    <PROTOCOLS>
        <Add mover="pack" />
        <Add mover="min_cart" />
    </PROTOCOLS>
```

	$> cp inputs/packing.xml .
	$> rosetta_scripts.default.linuxgccrelease -s 1ubq.pdb -parser:protocol packing.xml -out:prefix packing_ -nstruct 5

In the tracer output you should now see that both the PackRotamersMover and MinMover are running. We added the -nstruct 5 to produce five output structures. The Rosetta packer is stochastic, so different runs through the protocol should result in slightly different results. However, for repacking only (as opposed to design) the packer is rather convergent, so most structures should find about the same final conformation.

If you look at the scores of the packing run in comparison to the minimize run, you should see that the extra packing step allows us to sample a lower energy structure (about -190 REU versus -155 REU). Looking at the structures, you'll notice that they're mostly the same - especially in the core of the protein - but some of the surface sidechains have moved much more than they have from minimization only. 


Conclusion
----------

This tutorial was intended to give you a brief introduction in creating an XML protocol. The process we went through is similar to how most RosettaScripts developers write an XML file from scratch: Build up a protocol iteratively, starting with a simple protocol and progressively adding different and more complex stages. For each stage, have an idea about the effect you wish to accomplish, and then scan the documentation for existing movers/filters/task operations/etc. which will accomplish it. This may involve multiple RosettaScripts objects, due to movers which need as parameters other movers which need filters which need task operations (which need ...)

There are, of course, many more RosettaScripts objects than we have discussed, most of which should be covered in the RosettaScripts documentation. There are also additional sections of the XML, which are used for more specialized applications. (For example, ligand docking.) 

A final note - even if you can create an XML from scratch, it may be easier not to. If you already have an example XML that does something close to what you want to do, it's probably easier to start with that XML, and alter it to add in the functionality you want.

The hard part is not necessarily in putting together the XML, but in determining the optimal protocol (the logical steps) you should use to accomplish your modeling goals, and then in benchmarking that protocol to make sure it does what you hoped.

Troubleshooting
---------------

RosettaScripts is sensitive to mis-matched tags. If you forget to close a tag, or omit the ending slash on what is supposed to be a standalone tag, RosettaScripts will exit with a (possibly uninformative) error message. If you get something like Error: Tag::read - parse error", this means there is a syntax error in your XML. The recommended way of debugging it is to make a copy of the script, and progressively portions of the XML file until you get a script that works. (Or at least is able to be parsed.) It is then likely that the source of the error is in the portion of the XML which you deleted.
