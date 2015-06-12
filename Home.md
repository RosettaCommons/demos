Demos
=====

These demos are designed to guide users through sample procedures in computational modeling from the point of view of solving a specific problem.
Each one includes:

1. An introduction to the task at hand.
2. Detailed step-by-step instructions on how to run the demo.
3. All the input data needed to run the demo.
4. Scripts containing the exact commands needed to run the demo.

Downloading the demos
=====================

All the demos documented here are included in the `demos` directory at the root of the weekly releases of Rosetta.
Weekly releases are available to academic users free of charge.
If you are a [[Rosetta Developer]] (i.e. you have signed the RosettaCommons developer's agreement and have access to the Rosetta source code repositories on GitHub), you can download all the demos by checking out the `demos.git` repository:

    git clone git@github.com:RosettaCommons/demos.git

Demos
=====

How-Tos
-------

How-to demos demonstrate the best way to solve a particular problem using the Rosetta framework.
Each demo defines a particular problem, and as new and improved ways to solve those problems are developed, the demos are updated by members of the Rosetta community.
These demos are located in the `demos/public` directory.

<<LinkDemos(public)>>

Protocol Captures
-----------------

Protocol captures are demos from published papers. 
They aren't meant to show the best way to solve problems in the current version of Rosetta, they are meant to show published solutions to problems that were addressed by members of the Rosetta community.
The purpose of these protocol captures is both to serve as a historical record and to assist those trying to reproduce past results.
These demos are located in the `demos/protocol_capture` directory.

<<LinkDemos(protocol_capture)>>

<!--- BEGIN_INTERNAL --->

Adding new demos
================

To add new demos, the first thing you need to do is become a [[Rosetta developer]].
Then you will be able to check out the `demos.git` repository:

    git clone git@github.com:RosettaCommons/demos.git

The demos are organized into three directories:

* `demos/public`:

  For demos that are meant to show the best way to solve a particular problem.
These demos may be updated by the community as new ways to solve these problems are developed.

* `demos/protocol_capture`:

  For demos that are associated with published papers and tht demostrate the specific algorithm described in that paper.
These demos are static, may only work with previous old of rosetta, and meant to serve more as a historical records.

* `demos/pilot`:

  For demos that aren't meant to be included in the weekly releases yet.  I can't think of any good reasons for putting new demos here.

Each demo should go in its own directory within one of these three directories.
So to add a new demo, the first step is to create a descriptively named directory for it in the proper location.
For example, this is how you'd make a public demo called `my_demo`:

    cd demos/public
    mkdir my_demo

The only file that absolutely has to be in this directory is `README.md`.
A link to this file will automatically be added to the list of demos on this page.
Your readme should contain your name and email address, any relevant citations, a description of the problem your demo solves, links to all the scripts and input data your demo uses, and step-by-step instructions on how to run your demo.

Scripts and input files should of course be included in the directory as well.
How you organize these is up to you.
If you don't have many files, maybe just put everything in one directory.
If you have lots of files, maybe organize them into subdirectories.
Whatever makes the most sense for your demo.

Once you have finished writing your demo and have made sure that it runs properly, commit your changes and push them like usual: `git commit && git push`.
At this point, you can now edit your readme using the online Gollum interface available at the [[internal documentation site]].
However, only your readme can be maintained in this way.
You have to use git if you want to make any changes to your demo scripts or input data.

Under Development
-----------------

If you want to prevent a demo from being published to the static wiki until you've finished developing it or written a paper on it or something, put it in `demos/pilot`.
It will be visible here, but not on the public site.

<<LinkDemos(pilot)>>

<!--- END_INTERNAL --->
