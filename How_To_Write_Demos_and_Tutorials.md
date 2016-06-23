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
NOTE: If you add a protocol capture to the protocol_capture directory, you should add a copy to the public/ directory,
to serve as a version which can be updated to reflect best practices as Rosetta changes. 

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

In particular, it should have at least one line detailing the command that is to be executed to run your demo, starting with a dollar sign and greater than sign, which will be used to create the demo test for your demo. This will be parsed out of your demo's .md file and executed with the following substitutions:



Also, your flags files should be accompanied by a flags file of the same name with a ".short" suffix; this flags file will be used for the test preferentially over your actual flags file, and so may be used to minimize runtime at the expense of comprehensive demonstrativity.

 Using greater than signs (>) alone will not cause these commands to run, and can be used to indicate optional or supplementary commands to users of your demo. For example,

 > echo "This command will not be tested"
$> echo "This command will be tested."

Scripts and input files should of course be included in the directory as well; anything that Rosetta cannot generate should go into a folder specific to the inputs to the Rosetta job.
How you organize these is up to you.
If you don't have many files, maybe just put everything in one directory.
If you have lots of files, maybe organize them into subdirectories.
Whatever makes the most sense for your demo.
However, only .md files at the root of your demo directory will be tested; include cd commands to have the test navigate around to any subdirectories you may have..

###Keywords
Demos must be keyworded; that is, every MD file in your demo needs a line that reads something like
KEYWORDS: UTILITIES GENERAL
with at least two keywords drawn from the keywords.txt located in the demos directory. Omitting this line, omitting the keywords, or using unapproved keywords will cause the demo test to error out.


Once you have finished writing your demo and have made sure that it runs properly, commit your changes and push them like usual: `git commit && git push`.
At this point, you can now edit your readme using the online Gollum interface available at the [[internal documentation site]].
However, only your readme can be maintained in this way.
You have to use git if you want to make any changes to your demo scripts or input data.

A few days after you push you demo to the `demos` repository, your demo will 
become available from this website.  A link to it will automatically be added 
to the home page under the section indicated by your keywords.
 To make your demo easier to find, spend a few minutes browsing the 
[[documentation wiki|https://www.rosettacommons.org/docs/latest]] and adding 
links to any relevant pages you find.  The application section in particular 
would benefit from having lots of links to demos.  Note that (for technical 
reasons — Gollum gets really slow when there are too many pages in the wiki) 
the demos wiki is actually a whole different website than the documentation 
wiki.  So you have to use external links to link between the two wikis.
