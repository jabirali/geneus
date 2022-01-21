# GENEUS: General Non-Equilibrium Usadel Solver

You are currently viewing an EXPERIMENTAL branch, which will eventually
become GENEUS v2. This branch is not recommended for production use.
Please see the `master` branch for the stable version (GENEUS v1).

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
### Platforms and compilers
In contrast to GENEUS v1, this branch only supports a single `cmake` build
type, a single platform (Ubuntu x86-64), and a single compiler (GFortran).
Since the initial public release in 2018, I haven't had much time to maintain
the code (I don't use it that much anymore), and different bugs have surfaced
on different compilers and systems. By narrowing down the debugging to one
compilation target, the future maintenance of GENEUS becomes much easier.

### Speed versus robustness
One drawback of the above is performance. In my experience, GFortran can
produce 3-4x faster code with more aggressive optimization options (as
used by GENEUS v1), and Intel Fortran could produce 3-4x faster code than
that again, which in total could yield a speedup of an order of magnitude.

Thus, GENEUS v2 can be expected to be much slower than GENEUS v1. On the
other hand, less aggressive optimization also means that I don't have to
update the code every time a new compiler version is released, and that
I don't have to deal with different bugs across different compile targets.
This makes the code more robust and reliable – i.e. it should "just work".

### Code style
I've attempted to make the code conform more to the conventions used
in modern Fortran code, e.g. standard file extensions and line lengths.
These updates were performed with the help of the `fprettify` command.
Hopefully, this should make the code easier to read and modify by others.
