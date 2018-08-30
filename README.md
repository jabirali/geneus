This is a set of numerical programs for solving the Usadel diffusion equation in
superconducting nanostructures, both in and out of equilibrium. The programs are
quite flexible, being able to treat systems with e.g. magnetic elements, spin-orbit
coupling, spin-flip scattering, spin-orbit scattering, strongly polarized magnetic
interfaces, voltage gradients, and temperature gradients. The suite also contains
specialized programs for calculating the critical temperature and phase diagrams
of all these hybrid structures. The programs are also built to be user-friendly:
they are configured using a simple plaintext configuration file, which can include
mathematical expressions to initialize the physical properties of the system, and 
the output data is saved as plaintext files that can be directly imported in 
programs such as Gnuplot and Matlab. Finally, the code itself is written in 
modern object-oriented Fortran (2008+), making it both simple and efficient.

This software was developed by [Jabir Ali Ouassou](https://github.com/jabirali)
during his doctoral research in theoretical physics. The research was supervised
by [Prof. Jacob Linder](https://folk.ntnu.no/jacobrun/) at the 
[Center for Quantum Spintronics](https://www.ntnu.edu/quspin), NTNU, Norway.
The software itself is released under an MIT open-source licence; this basically means that
you are free to use it for any purpose, as long as you give credit where appropriate.

For more information, see the [documentation](doc/html/page/index.html).
