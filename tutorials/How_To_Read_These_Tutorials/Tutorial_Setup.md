#How To Read These Tutorials

KEYWORDS: CORE_CONCEPTS GENERAL

Written by Frank Teets
Last Modified Jun 21 2016

These tutorials were written such that a completely new user reading them [[in order|Home]] should gain a working understanding of the basic mechanics of Rosetta. There are many more protocols in Rosetta besides the basic ones the introductory tutorials cover. Many of these have [[demos available|demos-by-category]]. 

##Before Running Any of the Other Tutorials

To run these tutorials, you should have Rosetta installed and compiled. The [[Installation and Building|install_build]] tutorial 
should take you through the process, if it is not already installed. Please verify that the `<path_to_Rosetta_directory>/main/bin` directory contains executables appropriate to your installation. 

The tutorials and demos provided with Rosetta are written to assume that the following environment variables are set:

> env $ROSETTA3=the path to your Rosetta/source directory

> env $ROSETTA3_DB=the path to your Rosetta/database directory

> env $ROSETTA_TOOLS=the path to your Rosetta/tools directory

##Do The Following for Each Tutorial

In order for the hands-on portions of these tutorials to function correctly, you must make your current working directory that of the tutorial you want to run; i.e. for this tutorial, your current working directory must be (i.e. "you must be in") `${ROSETTA3}/demos/tutorials/How_To_Read_These_Tutorials`. Therefore, `cd` into the directory for the tutorial you want to run; see the [[Working with Rosetta|working_with_rosetta]] tutorial for a brief description of how to use `cd`.

You must also change the names of the executables to reflect the version that is compiled on your system. When in doubt, remove the final ".(os)(compiler)(mode)" element from the name of the executable and tab-complete to see what version you have compiled, then use that suffix in all tutorials.

If you have previously run a tutorial and wish to re-run it, make sure you delete any outputs created by the tutorials (e.g. the contents of the tutorial_output folder). By default, Rosetta will abandon a job rather than overwrite an existing output file, so the output directory must be clean for the tutorials to properly execute. One recommendation is to make a copy of the tutorial directory and run the examples there, which will make resetting the tutorial directory to a clean state easier.

