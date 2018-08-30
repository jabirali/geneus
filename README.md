# Superconducting spintronics simulation software

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
during his doctoral research. The research was supervised
by [Prof. Jacob Linder](https://folk.ntnu.no/jacobrun/) at the 
[Center for Quantum Spintronics](https://www.ntnu.edu/quspin), NTNU, Norway.
The software is released under the 
[MIT open-source licence](https://github.com/jabirali/DoctorCode/blob/master/LICENSE.md);
this basically means that you are free to use it for any purpose, 
as long as you give credit where appropriate.

The project also relies on the following Fortran libraries:

 * [`bvp_solver`](http://cs.stmarys.ca/~muir/BVP_SOLVER_Webpage.shtml),
   which is used to numerically solve boundary value problems;
 * [`fparser`](http://fparser.sourceforge.net/),
   which is used to parse mathematical expressions in configuration files;
 * [`pchip`](https://people.sc.fsu.edu/~jburkardt/f_src/pchip/pchip.html),
   which is used for one-dimensional interpolation.

The first is available under a mixture open-source licences, as indicated
in the source files. The latter two are available under the BSD
and LGPL licences, respectively. Their sources are located under `src/external/`.

For more information, see the 
[documentation](https://jabirali.github.io/DoctorCode/html/page/index.html).
