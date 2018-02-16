!> Author:   Jabir Ali Ouassou
!> Category: Programs
!>
!> This program calculates the phase-diagram of a one-dimensional superconducting structure.

program phase_p
  use :: structure_m
  use :: stdio_m
  use :: math_m

  !--------------------------------------------------------------------------------!
  !                                GLOBAL VARIABLES                                !
  !--------------------------------------------------------------------------------!

  ! Declare the superconducting structure
  type(structure)                 :: stack

  ! Declare program control parameters
  integer,      parameter         :: bootstraps = 10
  integer,      parameter         :: iterations = 10
  real(wp),     parameter         :: threshold  = 1e-8_wp

  ! Declare variables used by the program
  real(wp)                        :: flow
  real(wp)                        :: init



  !--------------------------------------------------------------------------------!
  !                           INITIALIZATION PROCEDURE                             !
  !--------------------------------------------------------------------------------!

  ! Redefine stdout and stderr 
  stdout = output('output.log')
  stderr = output('error.log')

  ! Construct the material stack
  stack = structure()

  ! Use the fixpoint method for selfconsistency iterations
  call stack % cmap('selfconsistency', 1)

  ! Find out what gap the user has initialized the system to
  init = stack % gap()

  ! Reset the states of the propagators throughout the stack
  call stack % initialize

  ! Bootstrap the material states while locking the gap
  call stack % converge(threshold = threshold, iterations = bootstraps, bootstrap = .true.)



  !--------------------------------------------------------------------------------!
  !                          PHASE DIAGRAM EVALUATION                              !
  !--------------------------------------------------------------------------------!

  ! Update the material states
  call stack % converge(iterations = iterations, prehook = prehook)

  ! Calculate the gap changes
  flow = stack % gap() / init

  ! Write out the final results
  call finalize



  !--------------------------------------------------------------------------------!
  !                                 SUBROUTINES                                    !
  !--------------------------------------------------------------------------------!

contains
  impure subroutine prehook
    ! Write out status information.
    flow = stack % gap() / init
    call status_body('Gap flow', flow)
  end subroutine

  impure subroutine finalize
    ! Status information
    call status_head('PHASE DIAGRAM')
    call status_body('Gap flow', flow)
    call status_foot

    ! Write the result to file
    call dump('flow.dat', flow)

    ! Close output files
    close(stdout)
    close(stderr)
  end subroutine
end program
