!> Author:   Jabir Ali Ouassou
!> Category: Programs
!>
!> This program calculates the critical temperature of an arbitrary superconducting hybrid structure
!> for a given set of physical parameters. In order to obtain the critical temperature as a function
!> of these parameters, the program has to be invoked multiple times with different input parameters.
!> The structure is constructed based on a configuration file which should be provided as the first 
!> command-line argument, and the results are then written to the output file 'critical.dat'. 

program main
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
  integer,      parameter         :: bootstraps = 20
  integer,      parameter         :: iterations = 10
  real(wp),     parameter         :: threshold  = 1e-08_wp
  real(wp),     parameter         :: initgap    = 1e-06_wp

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

  ! Disable all selfconsistency iterations
  call stack % cmap('selfconsistency', 0)

  ! Initialize the stack to a barely superconducting state
  call stack % initialize(cx(initgap))

  ! Bootstrap the material states at zero temperature
  call stack % converge(threshold = threshold, iterations = bootstraps, bootstrap = .true.)

  ! Save the current state of the materials
  call stack % save



  !--------------------------------------------------------------------------------!
  !                           BINARY SEARCH PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Use the fixpoint method for selfconsistency iterations
  call stack % cmap('selfconsistency', 1)

  do n = 1,bisections
    ! Set the temperature of the materials
    call stack % cmap('temperature', critical)

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
    ! Write out final results.

    ! Status information
    call status_head('CRITICAL TEMPERATURE')
    call status_body('Result', critical)
    call status_foot

    ! Write the critical temperature to file
    call dump('critical.dat', critical)

    ! Close output files
    close(stdout)
    close(stderr)
  end subroutine
end program
