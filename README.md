# GENEUS: General Non-Equilibrium Usadel Solver

You are currently viewing an EXPERIMENTAL branch, which is not recommended for
production use. Please checkout the `master` branch for the stable version.

In contrast to GENEUS v1, this branch only supports a single `cmake` build
type, a single platform (Ubuntu x86-64), and a single compiler (GFortran).
Moreover, this branch is quite a bit slower than GENEUS v1; it disables
most aggressive compilation options (e.g. link-time optimization, unsafe
math optimizations, and stack arrays), while enabling more reliability
(less work to port to new GFortran versions, clearer error messages).

To compile the project:

	./configure
	cd build
	make -j

The executables will then end up in the project subfolder `bin`. See
the documentation for more information on how to *use* these binaries.
