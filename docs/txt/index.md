title:  Overview
author: Jabir Ali Ouassou
date:   2018-08-30

@NOTE The documentation is incomplete, and parts of it are already out-of-date.
      The plan is to amend it with a discussion of the configuration file format as well as some example programs as soon as possible.

This is a set of numerical programs for solving the Usadel diffusion equation in superconducting nanostructures, both in and out of equilibrium.
The programs are quite flexible, being able to treat systems with e.g. magnetic elements, spin-orbit coupling, spin-flip scattering, spin-orbit scattering, strongly polarized magnetic interfaces, voltage gradients, and temperature gradients.
The suite also contains specialized programs for calculating the critical temperature and phase diagrams of all these hybrid structures.
The programs are also built to be user-friendly:
they are configured using a simple plaintext configuration file, which can include mathematical expressions to initialize the physical properties of the system, and the output data can be directly imported in e.g. Gnuplot or Matlab.
Finally, the code is written in modern object-oriented Fortran, making it both simple and efficient.

This software was developed by [Jabir Ali Ouassou](https://github.com/jabirali) during his doctoral research.
The research was supervised by [Prof. Jacob Linder](https://folk.ntnu.no/jacobrun/) at the [Center for Quantum Spintronics](https://www.ntnu.edu/quspin), NTNU, Norway.
The software itself is released under the [MIT open-source licence](https://github.com/jabirali/DoctorCode/blob/master/LICENSE.md);
this basically means that you are free to use it for any purpose, as long as you give credit where appropriate.
The source is available on [Github](https://github.com/jabirali/DoctorCode).

The project also relies on the external libraries [`bvp_solver`](http://cs.stmarys.ca/~muir/BVP_SOLVER_Webpage.shtml), [`fparser`](http://fparser.sourceforge.net/), and [`pchip`](https://people.sc.fsu.edu/~jburkardt/f_src/pchip/pchip.html).
The first is available under a mixture of open-source licences, as indicated in the source files themselves.
For convenience, these libraries have been bundled with this software page, and their sources are located under [src/external](https://github.com/jabirali/DoctorCode/tree/master/src/external).

You should start by checking if your system satisfies the [list of dependencies](01-dependencies.html), and then move on to the [compilation instructions](02-compilation.html).
After that, you can learn about the format of the configuration files used to describe relevant physical systems [here](03-config/index.html), and about each particular program in the suite [here](04-programs/index.html).

