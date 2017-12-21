!> Author:   Jabir Ali Ouassou
!> Category: Programs
!>
!> This program calculates the current-voltage relation of a one-dimensional superconducting structure.

program voltage_p
  use :: structure_m
  use :: stdio_m
  use :: math_m

  !--------------------------------------------------------------------------------!
  !                                GLOBAL VARIABLES                                !
  !--------------------------------------------------------------------------------!

  ! Declare the superconducting structure
  type(structure)                     :: stack
  type(conductor), target             :: na, nb

  ! Declare program control parameters
  real(wp), parameter                 :: threshold   = 1e-2
  real(wp), parameter                 :: tolerance   = 1e-8
  integer,  parameter                 :: iterations  = 1024

  ! Declare variables used by the program
  real(wp), dimension(:), allocatable :: voltage
  real(wp), dimension(:), allocatable :: current
  integer                             :: n



  !--------------------------------------------------------------------------------!
  !                           INITIALIZATION PROCEDURE                             !
  !--------------------------------------------------------------------------------!

  ! Redefine stdout and stderr 
  stdout = output('output.log')
  stderr = output('error.log')

  ! Construct the central material stack
  stack = structure()

  ! Ensure that the voltage has an effect
  call stack % cmap('equilibrium','F')

  ! Construct the surrounding conductors
  call na % construct()
  call nb % construct()

  ! Disable the conductors from updates
  call na % conf('order','0')
  call nb % conf('order','0')

  ! Thermalize the surrounding conductors
  na % temperature = stack % a % temperature
  nb % temperature = stack % b % temperature 

  ! Reinitialize the conductor states
  call na % initialize
  call nb % initialize

  ! Connect the conductors to the stack
  stack % a % material_a => na
  stack % b % material_b => nb

  ! Initialize the voltage and current
  allocate(voltage(iterations), current(iterations))
  call linspace(voltage, 0.0_wp, 2.00_wp)

  ! Start with a non-selfconsistent bootstrap procedure
  call stack % converge(threshold = threshold, bootstrap = .true.)



  !--------------------------------------------------------------------------------!
  !                           LINEAR SEARCH PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Calculate the charge current as a function of voltage difference
  current = 0
  do n=1,size(voltage)
    ! Update the voltage
    na % voltage = +voltage(n)/2
    nb % voltage = -voltage(n)/2

    ! Reset the states
    call na % initialize
    call nb % initialize

    ! Update the state
    call stack % converge(threshold = tolerance, prehook = prehook, posthook = posthook)

    ! Save the charge current to array
    current(n) = stack % a % supercurrent(0,1) + stack % a % lossycurrent(0,1)

    ! Write out the results
    call finalize
  end do



  !--------------------------------------------------------------------------------!
  !                                 SUBROUTINES                                    !
  !--------------------------------------------------------------------------------!

contains
  impure subroutine prehook
    ! Write out status information.
    call status_body('Voltage difference', voltage(n))
  end subroutine

  impure subroutine posthook
    ! Write results to output files.
    character(len=5) :: filename
    write(filename,'(f5.3)') voltage(n)
  end subroutine

  impure subroutine finalize
    ! Write the current-voltage relation to file
    call dump('current.dat', [voltage, current], ['Voltage           ', 'Charge current    '])
  end subroutine
end program
