!> Author:   Jabir Ali Ouassou
!> Category: Programs
!>
!> Calculates the current-voltage relation of a superconducting structure.

program voltage_p
    use :: structure_m
    use :: stdio_m
    use :: math_m

    !--------------------------------------------------------------------------!
    !                                GLOBAL VARIABLES                          !
    !--------------------------------------------------------------------------!

    ! Declare the superconducting structure
    type(structure)         :: stack
    type(conductor), target :: Ra, Rb

    ! Declare program control parameters
    real(wp), parameter :: threshold  = 1e-2
    real(wp), parameter :: tolerance  = 1e-5
    integer,  parameter :: iterations = 100

    ! Declare variables used by the program
    real(wp), dimension(:), allocatable :: voltage
    real(wp), dimension(:), allocatable :: current
    integer                             :: n

    !--------------------------------------------------------------------------!
    !                           INITIALIZATION PROCEDURE                       !
    !--------------------------------------------------------------------------!

    ! Construct the central material stack
    stack = structure()

    ! Ensure that the voltage has an effect
    call stack%cmap('nonequilibrium', 'T')

    ! Construct the surrounding conductors
    call Ra%construct()
    call Rb%construct()

    ! Disable the conductors from updates
    call Ra%conf('order', '0')
    call Rb%conf('order', '0')

    ! Thermalize the surrounding conductors
    Ra%temperature = stack%a%temperature
    Rb%temperature = stack%b%temperature

    ! Reinitialize the conductor states
    call Ra%initialize()
    call Rb%initialize()

    ! Connect the conductors to the stack
    stack%a%material_a => Ra
    stack%b%material_b => Rb

    ! Initialize the voltage and current
    allocate (voltage(iterations), current(iterations))
    call linspace(voltage, 1e-6_wp, 2.00_wp)

    ! Start with a bootstrap procedure (not self-consistent)
    call stack%converge(threshold=threshold, bootstrap=.true.)

    !--------------------------------------------------------------------------!
    !                           LINEAR SEARCH PROCEDURE                        !
    !--------------------------------------------------------------------------!

    ! Calculate the charge current as a function of voltage difference
    current = 0
    do n = 1, size(voltage)
        ! Update the voltage
        Ra%voltage = +voltage(n)
        Rb%voltage = -voltage(n)

        ! Reset the states
        call Ra%initialize()
        call Rb%initialize()

        ! Update the state
        call stack%converge(threshold=tolerance, prehook=prehook, posthook=posthook)

        ! Save the charge current to array
        current(n) = stack%a%supercurrent(0, 1) + stack%a%lossycurrent(0, 1)

        ! Write out the results
        call finalize()
    end do

    !--------------------------------------------------------------------------!
    !                                 SUBROUTINES                              !
    !--------------------------------------------------------------------------!

contains
    subroutine prehook
        call status_body('Voltage difference', voltage(n))
    end subroutine

    subroutine posthook
        character(len=5) :: filename
        write (filename, '(f5.3)') voltage(n)
    end subroutine

    subroutine finalize
        call dump('current.dat', [voltage, current], ['Voltage           ', 'Charge current    '])
    end subroutine
end program
