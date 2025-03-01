!> Author:   Jabir Ali Ouassou
!> Category: Programs
!>
!> Calculates general steady-state observables for a superconducting structure.

program converge_p
    use :: structure_m
    use :: stdio_m
    use :: math_m

    !--------------------------------------------------------------------------!
    !                                GLOBAL VARIABLES                          !
    !--------------------------------------------------------------------------!

    ! Declare the superconducting structure
    type(structure) :: stack

    ! Declare program control parameters
    real(wp), parameter :: threshold = 1e-2
    real(wp), parameter :: tolerance = 1e-8

    !--------------------------------------------------------------------------!
    !                           INITIALIZATION PROCEDURE                       !
    !--------------------------------------------------------------------------!

    ! Construct the superconducting structure
    stack = structure()

    ! Define output files
    stack%supercurrent  = output('supercurrent.dat')
    stack%lossycurrent  = output('lossycurrent.dat')
    stack%accumulation  = output('accumulation.dat')
    stack%distribution  = output('distribution.dat')
    stack%correlation   = output('correlation.dat')
    stack%magnetization = output('magnetization.dat')
    stack%density       = output('density.dat')

    !--------------------------------------------------------------------------!
    !                            CONVERGENCE PROCEDURE                         !
    !--------------------------------------------------------------------------!

    ! Bootstrap procedure (not self-consistent)
    call stack%converge(threshold=threshold, bootstrap=.true.)

    ! Convergence procedure (self-consistent)
    call stack%converge(threshold=tolerance)

    ! Write out the final results
    call finalize()

    !--------------------------------------------------------------------------!
    !                                 SUBROUTINES                              !
    !--------------------------------------------------------------------------!

contains
    subroutine finalize
        call status_head('STEADY STATE')
        call status_body('State difference', stack%difference())
        call status_body('Charge violation', stack%chargeviolation())
        call status_foot()
    end subroutine
end program
