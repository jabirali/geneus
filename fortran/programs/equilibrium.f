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

  ! Declare the superconducting structure
  type(structure) :: stack

  ! Declare program control parameters
  real(wp), parameter :: threshold = 1e-2
  real(wp), parameter :: tolerance = 1e-6

  ! Redefine stdout and stderr 
  stdout = output('output.log')
  stderr = output('error.log')

  ! Construct the superconducting structure
  stack = structure('materials.conf')

  ! Non-selfconsistent bootstrap procedure
  call stack % converge(threshold = threshold, bootstrap = .true.)

  ! Selfconsistent convergence procedure
  call stack % converge(threshold = tolerance, posthook = posthook)

  ! Write out the final results
  call finalize



  !--------------------------------------------------------------------------------!
  !                                 SUBROUTINES                                    !
  !--------------------------------------------------------------------------------!

contains
  impure subroutine posthook
    ! Write results to output files.
    use :: structure_m

    call stack % write_density('density.dat')
    call stack % write_current('current.dat')
    call stack % write_magnetization('magnetization.dat')
    call stack % write_gap('gap.dat')
  end subroutine

  impure subroutine finalize
    ! Write out the final results.
    use :: stdio_m

    ! Status information
    call status_head('EQUILIBRIUM')
    call status_body('State difference', stack % difference())
    call status_foot

    ! Close output files
    close(stdout)
    close(stderr)
  end subroutine
end program
