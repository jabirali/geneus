!> Author:   Jabir Ali Ouassou
!> Category: Programs
!>
!> This program calculates the critical temperature of a one-dimensional superconducting structure.

program critical_p
  use :: structure_m
  use :: stdio_m
  use :: math_m

  !--------------------------------------------------------------------------------!
  !                                GLOBAL VARIABLES                                !
  !--------------------------------------------------------------------------------!

  ! Declare the superconducting structure
  type(structure)                 :: stack

  ! Declare program control parameters
  integer,      parameter         :: bisections = 20
  integer,      parameter         :: bootstraps = 10
  integer,      parameter         :: iterations = 10
  real(wp),     parameter         :: threshold  = 1e-08_wp
  real(wp),     parameter         :: initgap    = 1e-02_wp

  ! Declare variables used by the program
  real(wp)                        :: minimum    = 0.00_wp
  real(wp)                        :: maximum    = 1.00_wp
  real(wp)                        :: critical   = 0.50_wp
  integer                         :: n          = 0



  !--------------------------------------------------------------------------------!
  !                           INITIALIZATION PROCEDURE                             !
  !--------------------------------------------------------------------------------!

  ! Redefine stdout and stderr 
  stdout = output('output.log')
  stderr = output('error.log')

  ! Construct the material stack
  stack = structure()

  ! Disable convergence acceleration
  call stack % cmap('boost', .false.)

  ! Initialize the stack to a barely superconducting state
  call stack % cmap('gap', initgap)

  ! Reset the states of the propagators throughout the stack
  call stack % initialize

  ! Bootstrap the material states while locking the gap
  call stack % converge(threshold = threshold, iterations = bootstraps, bootstrap = .true.)

  ! Save the current state of the materials
  call stack % save



  !--------------------------------------------------------------------------------!
  !                           BINARY SEARCH PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  do n = 1,bisections
    ! Set the temperature of the materials
    call stack % cmap('temperature', critical)

    ! Reinitialize at the new temperature
    call stack % initialize

    ! Load the saved material states
    call stack % load

    ! Update the material states
    call stack % converge(iterations = iterations, prehook = prehook)

    ! Update the critical temperature estimate
    if (stack % gap() >= initgap) then
      minimum = critical
    else
      maximum = critical
    end if
    critical = (minimum + maximum)/2
  end do

  ! Write out the final results
  call finalize



  !--------------------------------------------------------------------------------!
  !                                 SUBROUTINES                                    !
  !--------------------------------------------------------------------------------!

contains
  impure subroutine prehook
    ! Write out status information.
    call status_body('Temperature', critical)
    call status_body('Bisection',   n)
  end subroutine

  impure subroutine finalize
    ! Status information
    call status_head('CRITICAL TEMPERATURE')
    call status_body('Result', critical)
    call status_foot

    ! Write the result to file
    call dump('critical.dat', critical)

    ! Close output files
    close(stdout)
    close(stderr)
  end subroutine
end program
