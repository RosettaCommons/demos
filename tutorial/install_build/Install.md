# Downloading Rosetta

For any academic or commercial use, you need to [request a license](https://els.comotion.uw.edu/express_license_technologies/rosetta). Obtaining a license is free for academic users. After you obtained the license, you can [download Rosetta](https://www.rosettacommons.org/software/license-and-download). Make sure you download the version correspinding to the license you have. When you click, you can see the latest release as well as several weekly releases. The weekly release contains all the latest changes and but may lack documentations for certain new features. When you click on a version you can see that there are a number of options to download. The "__source__" is what you want to download. If you use any noncanonical amino acids, download the library of NCAAs as well.

# Installing Rosetta

The downloaded file is in form of [tar archive](https://en.wikipedia.org/wiki/Tar_(computing)) with .tgz extension. In a linux or mac, you can untar/uncompress the file by either double clicking on it or run this command in your terminal:
```
$tar -xvzf Rosetta[releasenumber].tar.gz
```
Unfortunatley,  currently there is no support for the whole Rosetta on Windows.  [Dual booting](https://en.wikipedia.org/wiki/Multi-booting) or [virtual machines](https://en.wikipedia.org/wiki/Virtual_machine) running Linux/MacOS are options. 

# Compiling Rosetta

Open the file that you unzipped and navigat through the folders: Rosetta -> main -> source or use the following bash command:
```
$cd Rosetta/main/source
```
You can see that currently the directory is empty. The reason is that the file you downloaded it comes in form of a collection of instructions, which is called [source code](https://en.wikipedia.org/wiki/Source_code). In order to be able to run or execute Rosetta, you need to first compile the code.

Rosetta uses [SCons](http://www.scons.org/) to assist compiling, so you need to first [download](http://scons.org/pages/download.html) and install it. You also need a C++ compiler to be installed on your computer (see [Install a complier](#Install-a-Compiler)).

Now you can build Rosetta using this general command line (make sure you are in the source folder)

```
$./scons.py -j <number_of_processors_to_use> mode=release bin
```
-j is indicating how many processors you want to use. This number depends on your computer. For example the command below uses 20 processors to build:

```
$./scons.py -j 20 mode=release bin
```
Expect a long time for the compilation to finish, several hours on one processor.

Now look at your source folder. There are several new folders in it, including bin. You can now run Rosetta!

### More Options for Building

As you noticed, the command you ran had a "mode" option and that you mentioned "bin" specifically. There are several other flags and modes that you can use to build certain parts or features.

- modes:
    - mode=release is the faster mode. It is tested and ready to go.
    - mode= debug (or not mentioning any mode) is slower and is mostly used for development purposes.

- specifying which parts to build
    - empty: if you don't provide any location, by diffult the libraries will be built.
    - bin: complete compliation
    - bin/rosetta_scripts.default.linuxgccrelease: only compiles the mentioned file(s).

- extras
    
    - extras=static builds static binaries. This can be useful for copying and running the apps on other systems.
    - extras=graphics mode enables OpenGL graphics for those apps that support it.
    - extras=mpi compiles rosetta in mpi [**M**assage **P**assing **I**nterface](https://computing.llnl.gov/tutorials/mpi/#What) format (for those executables that allow mpi run). for example the code below compiles only the rosetta_scripts in mpi format and release mode using 5 processors:

```
$./scons.py -j 5 mode=release bin/rosetta_scripts.mpi.linuxgccrelease extras=mpi
```
- Compiler specification
    - by default scons build rosetta using GCC compiler. However you can specify what compiler you want to use and what version using "__cxx__". For example the command below builds Rosetta completely in release mode using clang compiler version 4.5 using 10 processors:

```
$./scons.py -j 10 mode=release bin cxx=clang cxx_ver=4.5
```
# Building Rosetta using the Rosetta Xcode project (Mac)
 If you are interested in working with Rosetta code, you can Build Rosetta using the Rosetta Xcode project.You can use it to build, run, debug, browse, and edit the source code. You can find the isntructions on how to use Xcode to build Rosetta [here](https://www.rosettacommons.org/docs/latest/build_documentation/Build-Documentation).

# PyRosetta Download and Installation

[PyRosetta](http://www.pyrosetta.org/) is an interactive Python-based interface to Rosetta, allowing users to create custom molecular modeling algorithms with Rosetta sampling and scoring functions using Python scripting. PyRosetta was written for Python 2.6. You can follow instructions ond download and install [here](https://www.rosettacommons.org/docs/latest/scripting_documentation/PyRosetta/PyRosetta) and [here](http://www.pyrosetta.org/dow).

# Public Clusters with Rosetta

As part of the XSEDE initiative, the [TACC/Stampede](https://www.rosettacommons.org/docs/latest/build_documentation/TACC) cluster has Rosetta and PyRosetta centrally installed for authorized users.

# Install a Compiler

For Macs, install the XCode development packages. Even though you won't be compiling Rosetta through XCode, installing it will also install a compiler. (Clang, for recent versions of MacOS.)

For Linux, you will want to install the compiler package from your package management system. For Ubuntu and similar systems, the package "build-essential" installed with a command like sudo apt-get install build-essential will set your system up for compilation.

# Trouble shooting

Please check [Rosetta Documentation](https://www.rosettacommons.org/docs/latest/build_documentation/Build-Documentation) for details of the two most common errors that will come up during installation and how to deal with them.



