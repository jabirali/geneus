# GENEUS: General Non-Equilibrium Usadel Solver

You are currently viewing an EXPERIMENTAL branch, which will eventually
become GENEUS v2. This branch is not recommended for production use.
Please see the `master` branch for the stable version [GENEUS v1][v1].

## Quick start
You need a recent version of GNU/Linux with GFortran v8.x and CMake v3.x
installed. Note that the code currently segfaults when compiled with
GFortran v9.x or higher, and resolving this is still work in progress.
Meanwhile, I would recommend that you either compile the code using an
older GFortran version, or that you download precompiled GENEUS binaries.

If you use Ubuntu, the dependencies listed above can be installed via:

	sudo apt install cmake gfortran-8 gfortran-8-doc
	sudo update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-8 8


After cloning the repository and opening it in a terminal, you can then
compile the software via the following commands:

	./configure
	cd build
	make

The executables will then end up in the project subfolder `bin`. See
the [documentation][docs] for more information on how to use these.

## Differences from v1
### Platforms and compilers
In contrast to GENEUS v1, this branch only supports a single `cmake` build
type, a single platform (Ubuntu x86-64), and a single compiler (GFortran).
Since the initial public release in 2018, I haven't had much time to maintain
the code, and different bugs have surfaced on different compilers and systems.
Narrowing down the debugging to one compilation target simplifies maintenance.

### Speed versus robustness
One drawback of the above is performance. In my experience, GFortran
produces 3-4x faster code with more aggressive compilation options (as used
in GENEUS v1), and Intel Fortran can produce 3-4x faster code than that again.
Thus, GENEUS v2 can be up to an order of magnitude slower than v1 was.

On the other hand, less aggressive options also means that less rigorous
testing is required after each compiler update. Moreover, switching the compile
target from `native` to `x86-64` means that binaries can be shared between
computers instead of each person having to recompile the code from scratch.

### Code style
I've attempted to make the code conform more to the conventions used
in modern Fortran code, e.g. standard file extensions and line lengths.
These updates were performed with the help of [fprettify][fp]. Hopefully,
this should make the code easier to extend by others.

[v1]: https://github.com/jabirali/geneus/tree/v1.0
[docs]: https://jabirali.github.io/geneus/html/page/index.html
[fp]: https://github.com/pseewald/fprettify
