# GENEUS: General Non-Equilibrium Usadel Solver

## Quick start
You need a Linux/Unix system with GFortran v8.x and CMake v3.x.  Note that the code is susceptible to segfaults when compiled with GFortran v9.x or newer, and might not compile at all using older versions or other compilers.  I would therefore recommend that you install specifically GFortran 8.x if possible. If you use Ubuntu, the dependencies can be installed directly via:

	sudo apt install cmake gfortran-8 gfortran-8-doc
	sudo update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-8 8

On newer versions of Ubuntu, you might have to [download GFortran v8.x](https://gfortran.meteodat.ch/download/x86_64/) from the official sources. After cloning the repository and opening it in a terminal, you can then compile the software via the following:

	./configure
	cd build
	make

The executables will then end up in the project subfolder `bin`. See the [documentation][docs] for more information on how to use these binaries for physics simulations.

Before running any simulations, it's recommended that you add

	ulimit -s unlimited

to your shell config (e.g. `~/.bashrc`). If not, limitations in stack size may result in crashes if you attempt to multiple simulations in parallel.

## Differences from v1
In contrast to GENEUS v1, the current version only supports a single `cmake` build type, a single platform (Ubuntu x86-64), and a single compiler (GFortran).  Narrowing down the debugging to one compilation target simplifies maintenance.  The binaries are now also portable (fully static x86-64 binaries), meaning that
they can be transferred between computers without requiring recompilation.

One drawback of the above is performance. By dropping support for Intel Fortran, and switching to less aggressive compilation flags, GENEUS v2 runs up to an order of magnitude slower than v1 did. If this becomes a bottleneck, you are of course free to experiment with other compilation options; please don't be shy to reach out with potential bugfixes and optimizations in that case.

I've attempted to make the code conform more to the conventions used in modern Fortran code, e.g. standard file extensions and line lengths. These updates were performed with the help of [fprettify][fp]. Hopefully, this should make the code easier to extend by others if needed.

[v1]: https://github.com/jabirali/geneus/tree/v1.0
[docs]: https://jabirali.github.io/geneus/html/page/index.html
[fp]: https://github.com/pseewald/fprettify
