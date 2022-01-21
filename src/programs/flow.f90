!> Author:   Jabir Ali Ouassou
!> Category: Programs
!>
!> Calculates the phase diagram of a superconducting structure.

program flow_p
    use :: structure_m
    use :: stdio_m
    use :: math_m

    !--------------------------------------------------------------------------!
    !                                GLOBAL VARIABLES                          !
    !--------------------------------------------------------------------------!

    ! Declare the superconducting structure
    type(structure) :: stack

    ! Declare program control parameters
    integer,  parameter :: bootstraps = 10
    integer,  parameter :: iterations = 10
    real(wp), parameter :: threshold  = 1e-8_wp

    ! Declare variables used by the program
    real(wp) :: flow
    real(wp) :: init

    !--------------------------------------------------------------------------!
    !                           INITIALIZATION PROCEDURE                       !
    !--------------------------------------------------------------------------!

    ! Construct the material stack
    stack = structure()

    ! Disable convergence acceleration
    call stack%cmap('boost', .false.)

    ! Find out what gap the user has initialized the system to
    init = stack%gap()
    flow = 1.0

    ! Reset the states of the propagators throughout the stack
    call stack%initialize()

    ! Bootstrap the material states while locking the gap
    call stack%converge(threshold=threshold, iterations=bootstraps, bootstrap=.true.)

    !--------------------------------------------------------------------------!
    !                          PHASE DIAGRAM EVALUATION                        !
    !--------------------------------------------------------------------------!

    ! Update the material states
    call stack%converge(iterations=iterations, prehook=prehook)

    ! Calculate the gap changes
    flow = stack%gap()/init

    ! Write out the final results
    call finalize()

    !--------------------------------------------------------------------------!
    !                                 SUBROUTINES                              !
    !--------------------------------------------------------------------------!

contains
    subroutine prehook
        flow = stack%gap()/init
        call status_body('Gap flow', flow)
    end subroutine

    subroutine finalize
        call status_head('PHASE DIAGRAM')
        call status_body('Gap flow', flow)
        call status_foot()

        call dump('flow.dat', flow)
    end subroutine
end program
