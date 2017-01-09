!> Author:   Jabir Ali Ouassou
!> Date:     2016-03-08
!> Category: Programs
!>
!> This program calculates various physical observables for a superconducting thin-film structure
!> in equilibrium, including the charge and spin currents, density of states, and order parameter.
!> The heterostructure is constructed based on the configuration file 'materials.conf',  which is
!> expected to be in the runtime directory. The output is written to files in the same directory.

program equilibrium
  use :: structure_m
  use :: stdio_m
  use :: math_m

  ! Declare program control parameters
  real(wp), parameter :: threshold = 1e-2
  real(wp), parameter :: tolerance = 1e-6

  ! Create the superconducting structure
  type(structure) :: stack
  stack = structure('materials.conf')

  ! Non-selfconsistent bootstrap procedure
  call stack % converge(threshold = threshold, bootstrap = .true.)

  ! Selfconsistent convergence procedure
  call stack % converge(threshold = tolerance, output = .true.)

  ! Status information
  call status_head('EQUILIBRIUM')
  call status_body('State difference', stack % difference())
  call status_foot
end program
