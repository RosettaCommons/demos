## Symmetry within RosettaScripts

* * How to set up the scripts to handle and score symmetric structures*

Also refere to [symmetry user's guide](https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/symmetry#How-to-adopt-your-protocol-to-use-symmetry).

Manyp proteins have symmetric structures, i.e. they are replicates of one *primary* unit. Rosetta can handle symmetry using the symmetry code and have many protocols that are adopted to symmetry. Before running, make sure your protocol is symmetry-friendly (If not, see [How to adopt your protocol to adopt symmetry](https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/symmetry#How-to-adopt-your-protocol-to-use-symmetry)).

In order to run a symmetric protocol, there are certain things you need:

*Step 1:*

Generating symmetry definition file which is provided for you in this tutorial in the inputs directory. Please refer to [symmetry tutorial](../symmetry/symmetry.md) for a detailed description of symmetry files and how to work with them. Look at the symmetry_example.pdb in the inputs directory. You can see that it has a three_fold symmetry. That's what the C3.symm file provided in the inputs/ directory is telling Rosetta to set. 

*Step 2:*

Setting up the script for symmetry. We have provided a script for you in the scripts directory, named symmetry.xml. Let's take a look at the script.

You can see that on top in the `<SCOREFXNS>` part, we have added a new line:

```
<SCOREFXNS>
    <sfx_symm weights="talaris2014" symmetric=1 />
</SCOREFXNS>
```
This is required to tell Rosetta that you need to use symmetry score.

Now, let's go to the `<MOVERS>`. In the first line we added:

```
<SetupForSymmetry name=add_symm definition="C3.symm" />
```
This is probably the most important line! Here, you are telling Rosetta that it is going to be ran in the C3 mode. It will generate the fold tree so that the provided subunit will be replicated based on the corresponding symmetry definition. You can use `ExtractAsymmetricUnit` to do the reverse and extract the assymetric subunit from a symmetric pose.

In the next two lines we have:

```
<SymPackRotamersMover name="symm_pack" scorefxn=sfx_symm/>
<SymMinMover name="symm_min" scorefxn=sfx_symm bb=0 chi=1 jump=ALL  />
```
These are symmetry-adopted versions of PackRotamersMover and MinMover. jump should be set to all to refine the symmetric degrees of freedom. 


There are other symmetric-aware movers that can also be used. You can find them in the [Rosetta Documentation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/RosettaScripts). The C7.symm file is also in the input directory if you want to use.


**Running and Outputs**

Now let's run the scripts. You are going to use inputs/symm_test.pdb file as an example. This is just ubiquitin structure that has been cleaned and prepared. Based on the `PROTOCOLS` section, here is the order of running:

```
<Add mover=add_symm/>
<Add mover=symm_pack/>
<Add mover=symm_min/>
```
You can run the script using this command: ($ROSETTA3 is the path to your Rosetta/main/source)

```
$> $ROSETTA3/bin/rosettascripts.linuxgccrelease @symm.options
```

In your output, you can see that these lines are printed in the log file: (you can also find the whole tracer in outputs/symm_c3.log)

```
core.conformation.symmetry.util: =================== SYM FOLD TREE, jump notation: =symfixed= *indep* #symdof# jump[=follows] ========================
S1(229)
|----#j1#------>4:Sub1A(1-76)
|----=j4=---->S2(230)----j2=1----->80:Sub2A(77-152)
\----=j5=---->S3(231)----j3=1---->156:Sub3A(153-228)
```

This is were the symmetry is being defined and the jumps between subunits are set based on the symmetry definition file you provided.

The structure should look like something like this:


Now, based on the things you learnt from the [symmetry tutorial](../symmetry/symmetry.md), try to see if you can get other types of symmetric structures. An example of a C7 symmetric ubiquitin is provided for you in outputs/symm_test_c7.pdb:



