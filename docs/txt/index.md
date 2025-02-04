title:  Overview
author: Jabir Ali Ouassou
date:   2018-08-30

### GENEUS: General Non-Equilibrium Usadel Solver
GENEUS (pronounced "genius") is a set of numerical programs for solving the Usadel diffusion equation in one-dimensional superconducting nanostructures, both in and out of equilibrium.
The programs are quite flexible, being able to treat systems with e.g. magnetic elements, spin-orbit coupling, spin-dependent scattering, strongly polarized magnetic interfaces, voltage gradients, and temperature gradients.
The suite also contains specialized programs for calculating the critical temperature and phase diagrams of all these hybrid structures.
The programs are also built to be user-friendly: they are configured using simple configuration files, which can include mathematical expressions to initialize the physical system, and the output is easily imported in e.g. Gnuplot or Matlab.
Finally, the code is written in modern object-oriented Fortran, making it both simple and efficient.

This software was developed by [Jabir Ali Ouassou](https://github.com/jabirali) during his doctoral research.
The research was supervised by [Prof. Jacob Linder](https://folk.ntnu.no/jacobrun/) at the [Center for Quantum Spintronics](https://www.ntnu.edu/quspin), NTNU, Norway.
The software itself is released under the [MIT open-source licence](https://github.com/jabirali/GENEUS/blob/master/LICENSE.md);
this basically means that you are free to use it for any purpose, as long as you give credit where appropriate.
The source is available on [Github](https://github.com/jabirali/GENEUS).

The project also relies on the external libraries [`bvp_solver`](http://cs.stmarys.ca/~muir/BVP_SOLVER_Webpage.shtml), [`fparser`](http://fparser.sourceforge.net/), and [`pchip`](https://people.sc.fsu.edu/~jburkardt/f_src/pchip/pchip.html).
These are available under a mixture of open-source licences, as indicated in the source files themselves.
For convenience, these libraries have been bundled with this software page, and their sources are located under [src/external](https://github.com/jabirali/GENEUS/tree/master/src/external).

For instructions on how to install GENEUS, continue to the [next part](01-installation.html) of the manual.
After that, the next section will guide you through the use of the programs themselves.
