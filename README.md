# GENEUS: General Non-Equilibrium Usadel Solver

You are currently viewing an EXPERIMENTAL branch, which will eventually
become GENEUS v2. This branch is not recommended for production use.
Please see the `master` branch for the stable version [GENEUS v1][v1].

## Quick start
You need a recent version of GNU/Linux with GFortran and CMake installed.
If you're using e.g. Ubuntu, these can be installed via:

	sudo apt install gfortran cmake

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
the code (I don't use it that much anymore), and different bugs have surfaced
on different compilers and systems. By narrowing down the debugging to one
compilation target, the future maintenance of GENEUS becomes much easier.

### Speed versus robustness
One drawback of the above is performance. In my experience, GFortran
produces 3-4x faster code with more aggressive compilation options (as used
in GENEUS v1), and Intel Fortran produces 3-4x faster code than that again,
which yields a net speedup of an order of magnitude. Thus, GENEUS v2 can be
expected to be noticeably slower than GENEUS v1.

On the other hand, less aggressive options also means that I don't have 
to test and update the code rigorously every time a new compiler version is
released. Most things should "just work", making the code more reliable. I've
also changed the build target from `native` to `x86-64`, while retaining a
fully static compile system, making the binaries themselves more portable.

### Code style
I've attempted to make the code conform more to the conventions used
in modern Fortran code, e.g. standard file extensions and line lengths.
These updates were performed with the help of [fprettify][fp]. Hopefully,
this should make the code easier to extend by others.

[v1]: https://github.com/jabirali/geneus/tree/v1.0
[docs]: https://jabirali.github.io/geneus/html/page/index.html
[fp]: https://github.com/pseewald/fprettify
