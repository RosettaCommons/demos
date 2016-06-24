#How To Write Demos and Tutorials

KEYWORDS: UTILITIES GENERAL

Authors: Frank David Teets(teetsf@gmail.com), Jared Adolf-Bryfogle, Kale Kundert
Adding new demos
================

To add new demos, the first thing you need to do is become a [[Rosetta developer]].
Then you will be able to check out the `demos.git` repository:

    git clone git@github.com:RosettaCommons/demos.git

The demos are organized into four directories:

* `demos/tutorials`:
  Tutorials are intended to be general introductions to Rosetta that discuss basic Rosetta concepts or widely used protocols.
The material discussed in tutorials should be of interest to a majority of Rosetta users. 
(The expectation is that the naive Rosetta user would work through all the material in the Rosetta tutorials directory.)

* `demos/public`:  
  For demos that are meant to show the best way to solve a particular problem. 
The majority of protocols will have their demos here. 
These demos may be updated by the community as new ways to solve these problems are developed.

* `demos/protocol_capture`:  
  For demos that are associated with published papers and that demonstrate the specific algorithm described in that paper.
These demos are static, may only work with previous old of Rosetta, and meant to serve more as a historical records.
They WILL NOT be tested on the test server.

NOTE: If you add a protocol capture to the protocol_capture directory, you should add a copy to the public/ directory,
to serve as a version which can be updated to reflect best practices as Rosetta changes. See below for how to test this version of the protocol capture.

* `demos/pilot`:  
  For demos that aren't meant to be included in the weekly releases yet.
(Generally, though, use of git branches instead of pilot is recommended.)

Each demo should go in its own directory within one of these three directories.
So to add a new demo, the first step is to create a descriptively named directory for it in the proper location.
For example, this is how you'd make a public demo called `my_demo`:

    cd demos/public
    mkdir my_demo

The only file that absolutely has to be in this directory is `README.md`.
A link to this file will automatically be added to the list of demos on this page.
Your readme should contain your name and email address, any relevant citations, a description of the problem your demo solves, links to all the scripts and input data your demo uses, and step-by-step instructions on how to run your demo.

### Demo and Tutorial Testing
Demos and tutorials are now tested on a pass/fail basis on the test server with each merge of Rosetta main into master.  In order to get your demo or tutorial tested on the server, you can do either of these:

 * 1) Edit your .md file to add testing hooks that indicate to test that particular command (described below) (Recommended)
 
 * 2) Include a file called ```command``` in your demo or tutorial directory in the [exact] style of an integration script. Note that for log output of each command, all log redirection should start with 'log'.  I is recommended to number your logs so that they are not overridden.

Each demo/tutorial will start the test in the respective demo/tutorial just like integration tests. 

#### Testing from Documentation MD File
__The command you want to test should start with a dollar sign and greater than sign__, which will be used to create the demo test for your demo. This will be parsed out of your demo's .md file using __'$>'__ at the beginning of a line. Extra white space and markdown syntax will be striped from the line.

 Using greater than signs (>) alone will not cause these commands to run, and can be used to indicate optional or supplementary commands to users of your demo. For example,

 	> echo "This command will not be tested"
	$> echo "This command will be tested."

##### Substitutions
If you require certain paths, you may use substititions which will be replaced with the full path when the test is run through the integration test script
The following substitutions will occur:
  (add these to your shell profile or export them for easy testing of command):
	
	$ROSETTA3 -> Rosetta/main/source
	$ROSETTA3_DB -> Rosetta/main/database
	$ROSETTA_TOOLS -> Rosetta/tools
	$ROSETTA_DEMOS -> Rosetta/demos

##### Rosetta Executables
Please use app.default.linuxgccrelease as your path to the app.  Any extension with both the platform and compiler will be parsed.  The default -should- be present.  

This will be parsed correctly:

	my_app.default.linuxgccrelease @flags

This WILL NOT be parsed: 

	my_app.default.[platform][compiler]release


##### Shortening your run
__Your run should go no longer than ~15 minutes__.  In order to aid in this, your flags files can be accompanied by a flags file of the same name with a ".short" suffix that will be used to minimize runtime during testing   This [or these] flags files will added to the end of your command-line (enabling overrides from the short file).   

This short file should NOT be a copy of the flags file and should only be used to shorten the run.  This is so that two flags files will not have to be maintained.  Here is a simple example, consider an nstruct of 25 to demonstrate someting in your demo.  The short file can be a one-line file that overrides nstruct to 1.  The command in your md file would reference ```@flags```, while the test server would look for the .short file and if found will add it to your command: ```@flags @flags.short```
	

##### Etc. 

Scripts and input files should of course be included in the directory as well; anything that Rosetta cannot generate should go into a folder specific to the inputs to the Rosetta job.
How you organize these is up to you.
If you don't have many files, maybe just put everything in one directory.
If you have lots of files, maybe organize them into subdirectories.
Whatever makes the most sense for your demo.
__However, only .md files at the root of YOUR demo directory will be tested__; include cd commands to have the test navigate around to any subdirectories you may have if you wish..

__Please only use 1 command per line in your MD file or command for testing!__

###Keywords and Organization
Demos __MUST__ be keyworded.
That is, each tutorial and demo MD should have keyword list as follows:

	KEYWORDS: STRUCTURE_PREDICTION LOOPS GENERAL DESIGN MEMBRANES

Keywords one and two denote organizational levels.  If an item does not have a level2, use GENERAL.  Any extra keywords are used for keyword search.
The approved keywords can be found in keywords.txt in the demo root directory.  Keywords are present in order to have one approved name for something, such as SYMMETRY vs SYMMETRIC or DOCK vs DOCKING.  Omitting this line, omitting the keywords, or using unapproved keywords will cause the demo test to error out.

###Documentation Wiki
Once you have finished writing your demo and have made sure that it runs properly, commit your changes and push them like usual: `git commit && git push`.
At this point, you can now edit your readme using the online Gollum interface available at the [[internal documentation site]].
However, only your readme can be maintained in this way.
You have to use git if you want to make any changes to your demo scripts or input data.

A few days after you push your demo to the `demos` repository, your demo will 
become available from this website.  A link to it will automatically be added 
to the home page under the section indicated by your keywords.
 To make your demo easier to find, spend a few minutes browsing the 
[[documentation wiki|https://www.rosettacommons.org/docs/latest]] and adding 
links to any relevant pages you find.  The application section in particular 
would benefit from having lots of links to demos.  Note that (for technical 
reasons — Gollum gets really slow when there are too many pages in the wiki) 
the demos wiki is actually a whole different website than the documentation 
wiki.  So you have to use external links to link between the two wikis.
