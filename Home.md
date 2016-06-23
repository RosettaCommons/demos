Demos, Tutorials, and Protocol Captures
=======================================

<!-- Manual TOC, as the auto one is too bulky -->
<p></p><div class="toc"><div class="toc-title">Table of Contents</div>
<ul><ul><li><a href="#obtaining-tutorial-materials">Obtaining tutorial materials</a></li>
<li><a href="#tutorials">Tutorials</a></li>
<li><a href="#demos">Demos</a></li>
<li><a href="#protocol-captures">Protocol Captures</a></li>
<!--- BEGIN_INTERNAL -->
<li><a href="#adding-new-demos">Adding new demos</a></li>
<li><a href="#under-development">Under Development</a></li>
<!--- END_INTERNAL -->
</ul></ul></div>
<br/>

#### [[Demos listed by keyword|tag-search]]

Obtaining tutorial materials
============================

The demos, tutorials, protocol captures, and all example inputs are provided with the full Rosetta distribution, under the demos/ directory. Rosetta is available for license (free of charge to academic users) at <https://www.rosettacommons.org/software>.

<!--- BEGIN_INTERNAL -->
For RosettaCommons users, the demos repository should be automatically downloaded by the get_rosetta.sh download script. Alternatively, RosettaCommons users can download the demos repository from GitHub. e.g.

    git clone git@github.com:RosettaCommons/demos.git 

<!--- END_INTERNAL -->

Tutorials
=========

These are introductory tutorials intended as a gentle introduction to Rosetta concepts, and using common functionality of Rosetta. For additional examples and information on using Rosetta, see the demos (below) or the [Rosetta documentation](https://www.rosettacommons.org/docs/latest/)

Full input files for the tutorials are located in the `demos/tutorials/` directory of the Rosetta distribution. 

<<LinkDemos(tutorials)>>

Demos
=====

These demos are designed to guide users through sample procedures in computational modeling from the point of view of solving a specific problem. 

Full input files for the demos are located in the `demos/public/` directory of the Rosetta distribution.

<<LinkDemos(public)>>

Protocol Captures
=================

Many papers using Rosetta are accompanied by a protocol capture - an example of how to use the protocol discussed in the paper. The protocol captures below aren't meant to show the best way to solve problems in the current version of Rosetta, instead they are meant to show published solutions to problems that were addressed by members of the Rosetta community. The purpose of these protocol captures is both to serve as a historical record and to assist those trying to reproduce past results. See the demos (above) for updated versions of most protocol captures.

Full input files for the protocol captures are located in the demos/protocol_capture/ directory of the Rosetta distribution.

<<LinkDemos(protocol_capture)>>

<!--- BEGIN_INTERNAL --->

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

Scripts and input files should of course be included in the directory as well.
How you organize these is up to you.
If you don't have many files, maybe just put everything in one directory.
If you have lots of files, maybe organize them into subdirectories.
Whatever makes the most sense for your demo.

Once you have finished writing your demo and have made sure that it runs properly, commit your changes and push them like usual: `git commit && git push`.
At this point, you can now edit your readme using the online Gollum interface available at the [[internal documentation site]].
However, only your readme can be maintained in this way.
You have to use git if you want to make any changes to your demo scripts or input data.

A few days after you push you demo to the `demos` repository, your demo will 
become available from this website.  A link to it will automatically be added 
to this page.  However, it will be hard to find (especially for new users) with 
just that.  To make your demo easier to find, spend a few minutes browsing the 
[[documentation wiki|https://www.rosettacommons.org/docs/wiki]] and adding 
links to any relevant pages you find.  The application section in particular 
would benefit from having lots of links to demos.  Note that (for technical 
reasons — Gollum gets really slow when there are too many pages in the wiki) 
the demos wiki is actually a whole different website than the documentation 
wiki.  So you have to use external links to link between the two wikis.

Under Development
=================

If you want to prevent a demo from being published to the static wiki until you've finished developing it or written a paper on it or something, put it in `demos/pilot`.
It will be visible here, but not on the public site.

<<LinkDemos(pilot)>>

<!--- END_INTERNAL --->
