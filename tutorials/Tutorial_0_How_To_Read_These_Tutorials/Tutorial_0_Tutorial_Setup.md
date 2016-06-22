#How To Read These Tutorials
KEYWORDS: CORE_CONCEPTS GENERAL
Written by Frank Teets
Last Modified Jun 21 2016
These tutorials were written such that a completely new user should be able to complete them in numeric order through Tutorial 10 and thereby gain an understanding of the basic mechanics of Rosetta. Prior to doing any of the later tutorials, run through Tutorials 1-10.
##Before Running Any of the Other Tutorials
Complete [Tutorial 1](https://github.com/RosettaCommons/demos/blob/XRW2016/tutorials/install_build.md) to install and compile rosetta; verify that the `<path_to_Rosetta_directory>/main/bin` directory contains executables appropriate to your installation.
##Do The Following for Each Tutorial
In order for the hands-on portions of these tutorials to function correctly, you must make your current working directory that of the tutorial you want to run; i.e. for this tutorial, your current working directory must be (i.e. "you must be in") `<path_to_Rosetta_directory>/demos/tutorials/tutorial_0`.
You must also change the names of the executables to reflect the version that is compiled on your system. When in doubt, remove the final ".(os)(compiler)(mode)" element from the name of the executable and tab-complete to see what version you have compiled, then use that suffix in all tutorials.
Therefore, `cd` into the directory for the tutorial you want to run; see [Tutorial 2](https://github.com/RosettaCommons/demos/blob/XRW2016/tutorials/Tutorial_2_Working_With_Rosetta/Tutorial_2.md) for a brief description of how to use `cd`.
If you have previously run a tutorial and wish to re-run it, make sure you delete the contents of the tutorial_output folder; by default, Rosetta will abandon a job rather than overwrite an existing output file, so the output directory must be clean for the tutorials to properly execute.


