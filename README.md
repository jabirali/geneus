# GENEUS: General Non-Equilibrium Usadel Solver

GENEUS is a set of numerical programs for solving the Usadel diffusion equation in one-dimensional superconducting nanostructures, both in and out of equilibrium.
The programs are quite flexible, being able to treat systems with e.g. magnetic elements, spin-orbit coupling, spin-dependent scattering, strongly polarized magnetic interfaces, voltage gradients, and temperature gradients.
The suite also contains specialized programs for calculating the critical temperature and phase diagrams of all these hybrid structures.
The programs are built to be user-friendly: they are configured using simple configuration files, which can include mathematical expressions to initialize the physical system, and the output is easily imported in e.g. Gnuplot or Matlab.
Finally, the code is written in modern object-oriented Fortran, making it both simple and efficient.

This software was developed by [Jabir Ali Ouassou](https://orcid.org/0000-0002-3725-0885) during his doctoral research.
The research was supervised by [Prof. Jacob Linder](https://folk.ntnu.no/jacobrun/) at the [Center for Quantum Spintronics](https://www.ntnu.edu/quspin), NTNU, Norway.
The software itself is released under the [MIT open-source licence](https://github.com/jabirali/geneus/blob/master/LICENSE.md);
this basically means that you are free to use it for any purpose, as long as you give credit where appropriate.

The project also relies on the external libraries [`bvp_solver`](http://cs.stmarys.ca/~muir/BVP_SOLVER_Webpage.shtml), [`fparser`](http://fparser.sourceforge.net/), and [`pchip`](https://people.sc.fsu.edu/~jburkardt/f_src/pchip/pchip.html).
These are available under a mixture of open-source licences, as indicated in the source files themselves.
For convenience, these libraries have been bundled with this software page, and their sources are located under [src/external](https://github.com/jabirali/geneus/tree/master/src/external).

For more information, see the 
[documentation](https://jabirali.github.io/geneus/html/page/index.html).
