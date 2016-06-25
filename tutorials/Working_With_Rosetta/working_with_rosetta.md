# Working With Rosetta

KEYWORDS: CORE_CONCEPTS GENERAL

Written by Frank Teets
Last Modified Jun 21 2016

[[_TOC_]]

##Required Software

To productively use Rosetta, a number of additional programs are need:

  * 1. A text editor capable of outputting plain text files; For cmd-line tools, most people use [Emacs](https://www.gnu.org/software/emacs/manual/html_node/emacs/index.html), [Vim](https://www.washington.edu/computing/unix/vi.html), [Nano](https://www.nano-editor.org/dist/v2.0/nano.html), or similar. For GUI-based tools, [Sublime Text](https://www.sublimetext.com) and [TextMate](https://macromates.com), are very good on mac and Gedit is included with Ubuntu Linux. Note that *word processors* like Microsoft Word and Libre/Open Office are unsatisfactory, as the additional formating they introduce will confuse Rosetta. 

 * 2. A molecular viewing tool. Rosetta does not include a way to actually visualize the output files it creates; a tool such as [PyMOL](https://www.pymol.org/) or [Chimera](https://www.cgl.ucsf.edu/chimera/) is necessary to examine the output PDBs. Rosetta does include a PyMOLObserver for directly viewing its output in PyMOL; instructions for setting up PyMOLObserver are found [here](https://www.rosettacommons.org/docs/latest/rosetta_basics/graphics-and-guis).

 * 3. A terminal. Unix and Apple users have suitable terminals installed by default. Windows users attempting to run Rosetta remotely on a Mac or Unix machine will require a tool like PuTTY to provide the necessary interface.

It may also help to familiarize yourself with the basics of [command line interfaces](https://bash.cyberciti.biz/guide/Main_Page). Particularly useful commands include:

	ls foo			displays a list of files in the directory foo
	cd foo			moves the current working directory to the directory foo
	mv foo bar 		moves all of foo into bar if bar is a directory or into a file called bar if not
	cp foo bar		moves all of foo into bar if bar is a directory or into a file called bar if not
	ln -s <path_to_bar> foo establishes a link called a symbolic link or symlink from foo to bar
	top 			lists the currently running processes, together with their process IDs. This is useful for how much computational resource your pro