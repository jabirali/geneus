!! This program calculates various equilibrium properties for a superconducting multilayer thin-film structure,
!! such as the density of states, charge currents, spin currents, and superconducting gap in the structure.
!! @TODO: This module is supposed to replace/supersede the older program 'density.f'.
!!
!! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
!! Created: 2016-03-08
!! Updated: 2016-03-08

program equilibrium
  use mod_structure

  type(structure) :: stack
  stack = structure('simulation.conf')

  do while (stack % difference() > 1e-4)
    call stack % update
    call stack % write_density('density.dat')
    call stack % write_current('current.dat')
    call stack % write_gap('gap.dat')
  end do
end program
