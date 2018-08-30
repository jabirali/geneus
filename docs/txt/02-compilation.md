title:  Compilation
author: Jabir Ali Ouassou
date:   2018-08-30



@NOTE
  These instructions should work for Unix-like operating systems with all the [recommended dependencies](01-dependencies.html) installed.
  This should include Linux/Unix/Mac systems, and possibly also Windows if you have e.g. [Cygwin](https://www.cygwin.com/) installed. 
  If you use other platforms or compilers, you may have to compile the project manually.

After downloading the source code, open a terminal, and go to the top directory of the project.
This should be the directory with the files `configure` and `CMakeLists.txt`.
To prepare for compilation, you should then run one of the following commands, depending on whether you wish to compile with GFortran or IFort:

    ./configure --gfortran
    ./configure --ifort

If you need to search for bugs, it can be quite useful to enable debug mode.
This essentially means that all compiler optimizations are disabled, and extra safety checks are performed at both compiletime and runtime.
As a consequence, the generated code will run several times slower, but will also provide better feedback if something goes wrong.
To enable debug mode, run one of these commands instead of the above:

    ./configure --gfortran-debug
    ./configure --ifort-debug

After configuration, change to the subdirectory `build`, and run `make` to compile:

    cd build
    make -j4

If the project compiles successfully, the executables should appear in the subdirectory `bin` of the top directory. 
You might want to add that directory to your system `PATH` variable, or move the executables to a system location like `/usr/local/bin/`, so that the programs can be run without specifying a full path.

This documentation has been automatically generated using [FORD](https://github.com/cmacmackin/ford).
If you modify the source code, and wish to update the documentation accordingly, run the following command:

    ./configure --docs

