title:  Equilibrium
author: Jabir Ali Ouassou
date:   2016-04-13

This program does something. Fill in details.

    ./bin/equilibrium

The boundary value problem solver used internally by the material objects works by starting with
an initial guess for what the solution should look like, and then tries to locally minimize the
implicit differential equations to find the proper solution. The better the initial guess is, the
better able the boundary value problem solver is to find a proper solution. By default, initial
states of all regular materials are set to BCS superconducting states, while the initial states
of halfmetallic ferromagnets are set to normal metallic states. However, for many multi-layer
structures, these are quite bad initial guesses. One workaround is to then start with a very
faint superconducting gap in the structure, and let it selfconsistently tend to the real one.

But for tough problems, you should use bootstrap binary instead.
