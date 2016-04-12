!! This program calculates various physical observables for a superconducting thin-film structure
!! in equilibrium, including the charge and spin currents, density of states, and order parameter.
!! The heterostructure is constructed based on the configuration file 'materials.conf',  which is
!! expected to be in the runtime directory. The output is written to files in the same directory.
!!
!! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
!! Created: 2016-03-08
!! Updated: 2016-03-23

program equilibrium
  use structure_m
  use stdio_m
  implicit none

  ! Create the superconducting structure
  type(structure) :: stack
  stack = structure('materials.conf')

  ! Main loop
  do
    ! Status information
    call status_head('UPDATING STATE')
    call status_body('State difference', stack % difference())
    call status_foot

    ! Update the material state
    call stack % update

    ! Write the results to files
    call stack % write_density('density.dat')
    call stack % write_current('current.dat')
    call stack % write_gap('gap.dat')

    ! Stop if we have convergence
    if (stack % difference() < 1e-6) then
      exit
    end if
  end do

  ! Status information
  call status_head('CONVERGED')
  call status_body('State difference', stack % difference())
  call status_foot
end program
