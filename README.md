# GENEUS: General Non-Equilibrium Usadel Solver

You are currently viewing an EXPERIMENTAL branch, which is not recommended for
production use. Please checkout the `master` branch for the stable version.

## Quick start
You need an up-to-date installation of GNU/Linux with GFortran and CMake
installed. If you're using e.g. Ubuntu, these tools can be installed via:

	sudo apt install gfortran cmake

After cloning the repository and opening it in a terminal, you can then
compile the software by typing the following commands into the terminal:

	./configure
	cd build
	make -j

The executables will then end up in the project subfolder `bin`. See
the documentation for more information on how to *use* these binaries.

## Differences from v1
### Platform support
In contrast to GENEUS v1, this branch only supports a single `cmake` build
type, a single platform (Ubuntu x86-64), and a single compiler (GFortran).
Since the initial public release in 2018, I haven't had much time to
maintain the code, and since different bugs tend to surface for each
target the code is compiled for, I've decided to narrow it down.

This also entails a different trade-off between speed and reliability.
The previous `Release` build was a factor 3x faster than the `Debug`
build when using GFortran, and with the Intel Fortran compiler the
result was typically 3x-4x faster than that. For GENEUS v2, I've kept
only one modified `Debug` build target, which can be up to an order
of magnitude slower than GENEUS v1. The trade-off here is reliability
and maintainability; disabling aggressive optimization options (e.g.
"unsafe" math optimizations, link-time optimizations, stack arrays)
drastically reduces the amount of bugs I see when compiling the code
with new compiler versions. Thus, the effort needed to keep the code
functional and reliable on modern systems is proportionally smaller.

### Code style
I've attempted to make the code conform more to the conventions used
in modern Fortran code, e.g. standard file extensions and line lengths.
These updates were performed with the help of the `fprettify` command.
Hopefully, this should make the code easier to modify by other authors.


