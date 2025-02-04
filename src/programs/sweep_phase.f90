!> Author:   Jabir Ali Ouassou
!> Category: Programs
!>
!> Calculates the current-phase relation of a superconducting structure.

program phase_p
    use :: structure_m
    use :: stdio_m
    use :: math_m

    !--------------------------------------------------------------------------!
    !                                GLOBAL VARIABLES                          !
    !--------------------------------------------------------------------------!

    ! Declare the superconducting structure
    type(structure)              :: stack
    type(superconductor), target :: Ra, Rb

    ! Declare program control parameters
    real(wp), parameter :: threshold  = 1e-2
    real(wp), parameter :: tolerance  = 1e-8
    integer,  parameter :: iterations = 51

    ! Declare variables used by the program
    real(wp), dimension(:), allocatable :: phase
    real(wp), dimension(:), allocatable :: current
    real(wp)                            :: critical
    integer                             :: n, m

    !--------------------------------------------------------------------------!
    !                           INITIALIZATION PROCEDURE                       !
    !--------------------------------------------------------------------------!

    ! Construct the central material stack
    stack = structure()

    ! Construct the surrounding superconductors
    call Ra%construct()
    call Rb%construct()

    ! DiRable the superconductors from updates
    call Ra%conf('order', '0')
    call Rb%conf('order', '0')

    ! Connect the superconductors to the stack
    stack%a%material_a => Ra
    stack%b%material_b => Rb

    ! Count the number of superconductors
    n = stack%superconductors()

    ! Depending on the number of enabled superconductors, we might get a 2πn
    ! periodicity instead of 2π periodicity, and need to account for this.
    m = (n + 1)*(iterations - 1) + 1
    allocate (phase(m), current(m))
    call linspace(phase, 1e-6_wp, 2*(n + 1) - 1e-6_wp)

    ! Start with a bootstrap procedure (not self-consistent)
    call stack%converge(threshold=threshold, bootstrap=.true.)

    !--------------------------------------------------------------------------!
    !                           LINEAR SEARCH PROCEDURE                        !
    !--------------------------------------------------------------------------!

    ! Calculate the charge current as a function of phase difference
    do n = 1, size(phase)
        ! Update the phase
        Ra%correlation = exp(((0.0, -0.5)*pi)*phase(n))
        Rb%correlation = exp(((0.0, +0.5)*pi)*phase(n))

        ! Reset the states
        call Ra%initialize()
        call Rb%initialize()

        ! Update the state
        call stack%converge(threshold=tolerance, prehook=prehook, posthook=posthook)

        ! Save the charge current to array
        current(n) = stack%a%supercurrent(0, 1)
    end do

    ! Calculate the critical current
    critical = maxval(abs(current))

    ! Write out the final results
    call finalize()

    !--------------------------------------------------------------------------!
    !                                 SUBROUTINES                              !
    !--------------------------------------------------------------------------!

contains
    subroutine prehook
        call status_body('Phase difference', phase(n))
    end subroutine

    subroutine posthook
        character(len=5) :: filename
        write (filename, '(f5.3)') phase(n)
    end subroutine

    subroutine finalize
        call status_head('CRITICAL CURRENT')
        call status_body('Result', critical)
        call status_foot()

        call dump('current.dat', [phase, current], ['Phase             ', 'Charge current    '])
        call dump('critical.dat', critical)
    end subroutine
end program
