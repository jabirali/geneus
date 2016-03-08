!! This program calculates various equilibrium properties for a superconducting multilayer thin-film structure,
!! such as the density of states, charge currents, spin currents, and superconducting gap in the structure.
!! @TODO: This module is supposed to replace/supersede the older program 'density.f'.
!!
!! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
!! Created: 2016-03-08
!! Updated: 2016-03-08

program equilibrium
  use mod_structure
  use mod_math

  type(structure) :: bilayer
  bilayer = structure('simulation.conf')

  do while (bilayer % difference() > 1e-4)
    call bilayer % update
    call bilayer % write_density('density.dat')
    call bilayer % write_current('current.dat')
    call bilayer % write_gap('gap.dat')
  end do
end program
