!> Author:   Jabir Ali Ouassou
!> Category: Programs
!>
!> This program calculates various physical observables for a superconducting thin-film structure
!> in equilibrium, including the charge and spin currents, density of states, and order parameter.
!> The heterostructure is constructed based on the configuration file 'materials.conf',  which is
!> expected to be in the runtime directory. The output is written to files in the same directory.

program main
  use :: structure_m
  use :: stdio_m
  use :: math_m

  ! Declare the superconducting structure
  type(structure) :: stack

  ! Declare program control parameters
  real(wp), parameter :: threshold = 1e-2
  real(wp), parameter :: tolerance = 1e-8

  ! Redefine stdout and stderr 
  stdout = output('output.log')
  stderr = output('error.log')

  ! Construct the superconducting structure
  stack = structure()

  ! Redefine output files
  stack % supercurrent  = output('supercurrent.dat')
  stack % lossycurrent  = output('lossycurrent.dat')
  stack % accumulation  = output('accumulation.dat')
  stack % correlation   = output('correlation.dat')
  stack % density       = output('density.dat')

  ! Non-selfconsistent bootstrap procedure
  call stack % converge(threshold = threshold, bootstrap = .true.)

  ! Selfconsistent convergence procedure
  call stack % converge(threshold = tolerance)

  ! Write out the final results
  call finalize



  !--------------------------------------------------------------------------------!
  !                                 SUBROUTINES                                    !
  !--------------------------------------------------------------------------------!

contains
  impure subroutine finalize
    ! Write out the final results.
    use :: stdio_m

    ! Status information
    call status_head('EQUILIBRIUM')
    call status_body('State difference', stack % difference())
    call status_body('Charge violation', stack % chargeviolation())
    call status_foot

    ! Close output files
    close(stdout)
    close(stderr)
  end subroutine
end program
