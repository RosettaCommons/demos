# Installing Rosetta

KEYWORDS: CORE_CONCEPTS GENERAL

*See also the Rosetta build/install documentation [here](https://www.rosettacommons.org/docs/latest/build_documentation/Build-Documentation).*

[[_TOC_]]

# Downloading Rosetta

For any academic or commercial use, you need to [request a license](https://els.comotion.uw.edu/express_license_technologies/rosetta). Obtaining a license is free for academic users. After you obtained the license, you can [download Rosetta](https://www.rosettacommons.org/software/license-and-download). Make sure you download the version corresponding to the license you have. When you click, you can see the latest numbered release as well as several weekly releases. Numbered releases (since Rosetta3.6) are simply weekly releases that have been specially labeled - all weekly releases pass the same suite of tests that the numbered releases do.

For many version, we offer both a source and a binary version. The binary version may allow you to skip the compilation stage, but are more limited in the platforms on which they work. The "__source__" distribution should be useful on all platforms on which Rosetta can run. (If you're interested in noncanonical amino acids, download the NCAA rotamer libraries as well.)

# Installing Rosetta

The downloaded file is in form of [tar archive](https://en.wikipedia.org/wiki/Tar_(computing)) with .tgz extension. In a linux or mac, you can untar/uncompress the file by either double clicking on it or run this command in your terminal:

```
> tar -xvzf rosetta[releasenumber].tar.gz
```

Unfortunately, currently there is no support for the whole Rosetta on Windows.  [Dual booting](https://en.wikipedia.org/wiki/Multi-booting) or [virtual machines](https://en.wikipedia.org/wiki/Virtual_machine) running Linux/MacOS are options. 

# Compiling Rosetta

Open the file that you unzipped and navigate through the folders: Rosetta -> main -> source or use the following bash command:

```
> cd rosetta*/main/source
```

If you downloaded the source bundle, you can see that the bin/ directory is currently empty. In order to be able to run Rosetta, you need to first compile the code.

To compile Rosetta you need a C++ compiler. Rosetta developers typically use [GCC](https://gcc.gnu.org/) or [Clang](http://clang.llvm.org/), although other standard-compliant compilers can be used. (See [Install a complier](#Install-a-Compiler) for more information on installing a compiler.)

Rosetta uses [SCons](http://www.scons.org/) as a build system. While Scons is available as a separate download, the Rosetta download includes a version, which is the recommended version to use in compiling Rosetta.

Now you can build Rosetta using this general command line (make sure you are in the source folder)

```
> ./scons.py -j <number_of_cores_to_use> mode=release bin
```

`-j` is indicating how many cores you want to use. This number depends on your computer. For example the command below uses 20 cores to build:

```
> ./scons.py -j 20 mode=release bin
```

Expect a long time for the compilation to finish, several hours on one core.

Now look at your source folder. There are several new folders in it, including bin. You can now run Rosetta!

### More Options for Building

As you noticed, the command you ran had a "mode" option and that you mentioned "bin" specifically. There are several other flags and modes that you can use to build certain parts or features.

- modes:
    - mode=release compiles with optimizations to produce a faster version of Rosetta.
    - mode=debug (or not mentioning any mode) includes additional checks which slows down Rosetta runs. It is mostly used for development and debugging purposes.

- specifying which parts to build
    - empty: if you don't provide any location, by default only the libraries will be built.
    - "bin": complete compilation of all applications in the bin/ directory 
    - "bin/rosetta_scripts.default.linuxgccrelease" or "rosetta_scripts": only compiles the mentioned application (multiple can be listed)

- extras
    - "extras=static": builds static binaries. This can be useful for copying and running the apps on other systems.
    - "extras=graphics": mode enables OpenGL graphics for those apps that support it.
    - "extras=opencl": enables GPU usage for those apps that support it
    - "extras=mpi" compiles Rosetta in MPI [**M**assage **P**assing **I**nterface](https://computing.llnl.gov/tutorials/mpi/#What) format (for those executables that support MPI runs). Running Rosetta in MPI mode requires (potentially non-trivial) edits to the site.settings files.

For example the code below compiles only the rosetta_scripts application in MPI format and release mode using 5 cores:

```
> ./scons.py -j 5 mode=release bin/rosetta_scripts.mpi.linuxgccrelease extras=mpi
```

**NOTE** when you build with different extras, the extension will change. For example if you use extras=mpi, you use rosetta_scripts.mpi.linuxgccrelease instead of rosetta_scripts.default.linuxgccrelease

- Compiler specification
    - by default scons builds Rosetta using the GCC compiler. However you can specify what compiler and version you wish to use by using "__cxx__". For example the command below builds Rosetta completely in release mode using clang compiler version 4.5 using 10 cores:

```
> ./scons.py -j 10 mode=release bin cxx=clang cxx_ver=4.5
```

# Building Rosetta using the Rosetta Xcode project (Mac)

If you are interested in working with Rosetta code, you can build Rosetta using the Rosetta Xcode project. You can use it to build, run, debug, browse, and edit the source code. You can find the instructions on how to use Xcode to build Rosetta [here](https://www.rosettacommons.org/docs/latest/build_documentation/Build-Documentation).

# PyRosetta Download and Installation

[PyRosetta](http://www.pyrosetta.org/) is an interactive Python-based interface to Rosetta, allowing users to create custom molecular modeling algorithms with Rosetta sampling and scoring functions using Python scripting. PyRosetta was written for Python 2.6. You can follow instructions to download and install PyRosetta [here](https://www.rosettacommons.org/docs/latest/scripting_documentation/PyRosetta/PyRosetta) and [here](http://www.pyrosetta.org/dow).

# Public Clusters with Rosetta

As part of the XSEDE initiative, the [TACC/Stampede](https://www.rosettacommons.org/docs/latest/build_documentation/TACC) cluster has Rosetta and PyRosetta centrally installed for authorized users.

# Install a Compiler

For Macs, install the XCode development packages. Even though you won't be compiling Rosetta through XCode, installing it will also install a compiler. (Clang, for recent versions of MacOS.)

For Linux, you will want to install the compiler package from your package management system. For Ubuntu and similar systems, the package "build-essential" can be installed with a command like `sudo apt-get install build-essential`.

# Trouble shooting

Please check the [Rosetta build documentation](https://www.rosettacommons.org/docs/latest/build_documentation/Build-Documentation) for details on common errors that will come up during installation and how to deal with them.



