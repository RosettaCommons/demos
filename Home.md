Demos
=====

These demos are designed to guide users through sample procedures in 
computational modeling from the point of view of solving a specific problem. 
Each one includes:

1. An introduction to the task at hand.
2. Detailed step-by-step instructions on how to run the demo.
3. All the input data needed by the demo.
4. Scripts that can be run as is to process the included input.

Downloading the demos
=====================

All the demos listed below should include links that allow you to download all 
the scripts and input data needed to run the demo.  If those links are missing, 
email the author of the demo and ask them to fix it.  Alternatively, if you 
are a [[Rosetta developer]] (i.e. you have signed the RosettaCommons developer's 
agreement and have access to the Rosetta source code repositories on GitHub), 
you can download all the demos by checking out the `documentation.git` 
repository:

    git clone git@github.com:RosettaCommons/documentation.git

Once you have done this, you should add public links (see the section on adding 
new demos below) then commit and push your changes.  

Adding new demos
================

To add new demos, the first thing you need to do is become a [[Rosetta 
developer]].  Then you will be able to checkout the `documentation.git` 
repository:

    git clone git@github.com:RosettaCommons/documentation.git

The demos are located in `application_documentation/demos`.  Create a new 
directory within this directory, and give it a descriptive name.  For example, 
this is how you'd make a demo called `my_demo`:

    cd documentation/application_documentation/demos
    mkdir my_demo

The only file that absolutely has to be in this directory is `readme.md`.  A 
link to this file will automatically be added to the list of demos below.  Your 
readme should contain your name and email address, any relevant citations, a 
description of the problem your demo solves, links to all the scripts and input 
data your demo uses, and step-by-step instructions on how to run your demo.

Scripts and input files should of course be included in the directory as well.  
How you organize these is up to you.  If you don't have many files, maybe just 
put everything in one directory.  If you have lots of files, maybe organize 
them into subdirectories.  Whatever makes the most sense for your demo.  
However, it is absolutely essential that your readme file includes Markdown 
links to all the scripts and input data that your demo uses.  If these links 
are missing, people who aren't Rosetta developers won't be able to access these 
files and won't be able to run your demo.

Once you have finished writing your demo and have made sure that it runs 
properly, commit your changes and push them like usual:

    git commit
    git push

At this point, you can now edit your readme using the online Gollum interface 
available at the [internal documentation site].  You can also continue to 
maintain your demo through the git interface if you prefer that.

Demos
=====

How-Tos
-------

How-to demos demonstrate the best way to solve a particular problem using the 
Rosetta framework.  Each demo defines a particular problem, and as new and 
improved ways to solve those problems are developed, the demos will be updated 
by members of the Rosetta community.

<<LinkDemos(public)>>

Protocol Captures
-----------------

Protocol captures are demos from published papers.  They aren't meant to 
show the best way to solve a certain problem in the current version of Rosetta, 
they are meant to show how the problem was solved by the authors of a 
publication at the time they addressed it.  The purpose of these protocol 
captures is both to serve as a historical record and to assist those trying to 
reporoduce past results.

<<LinkDemos(protocol_capture)>>

<!--- BEGIN_INTERNAL --->

Under Development
-----------------

If you want to prevent a demo from being published to the static wiki until 
you've finished developing it or written a paper on it or something, put it 
in `applications/demos/pilot`.  It will be visible here, but not on the public 
site.

<<LinkDemos(pilot)>>

<!--- END_INTERNAL --->
