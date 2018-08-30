!> Author:   Jabir Ali Ouassou
!> Category: Programs
!>
!> This program calculates the current-phase relation of a one-dimensional superconducting structure.

program phase_p
  use :: structure_m
  use :: stdio_m
  use :: math_m

  !--------------------------------------------------------------------------------!
  !                                GLOBAL VARIABLES                                !
  !--------------------------------------------------------------------------------!

  ! Declare the superconducting structure
  type(structure)                 :: stack
  type(superconductor), target    :: sa, sb

  ! Declare program control parameters
  real(wp), parameter             :: threshold   = 1e-2
  real(wp), parameter             :: tolerance   = 1e-8
  integer,  parameter             :: iterations  = 51

  ! Declare variables used by the program
  real(wp), dimension(:), allocatable :: phase
  real(wp), dimension(:), allocatable :: current
  real(wp)                            :: critical
  integer                             :: n, m



  !--------------------------------------------------------------------------------!
  !                           INITIALIZATION PROCEDURE                             !
  !--------------------------------------------------------------------------------!

  ! Redefine stdout and stderr 
  stdout = output('output.log')
  stderr = output('error.log')

  ! Construct the central material stack
  stack = structure()

  ! Construct the surrounding superconductors
  call sa % construct()
  call sb % construct()

  ! Disable the superconductors from updates
  call sa % conf('order','0')
  call sb % conf('order','0')

  ! Connect the superconductors to the stack
  stack % a % material_a => sa
  stack % b % material_b => sb

  ! Count the number of superconductors
  n = stack % superconductors()

  ! Depending on the number of enabled superconductors in the junction, we might get a
  ! 2πn-periodicity instead of a 2π-periodicity, and should therefore account for that.
  m = (n+1)*(iterations-1) + 1
  allocate(phase(m), current(m))
  call linspace(phase, 1e-6_wp, 2*(n+1) - 1e-6_wp)

  ! Start with a non-selfconsistent bootstrap procedure
  call stack % converge(threshold = threshold, bootstrap = .true.)



  !--------------------------------------------------------------------------------!
  !                           LINEAR SEARCH PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Calculate the charge current as a function of phase difference
  do n=1,size(phase)
    ! Update the phase
    sa % correlation = exp(((0.0,-0.5)*pi)*phase(n))
    sb % correlation = exp(((0.0,+0.5)*pi)*phase(n))

    ! Reset the states
    call sa % initialize
    call sb % initialize

    ! Update the state
    call stack % converge(threshold = tolerance, prehook = prehook, posthook = posthook)

    ! Save the charge current to array
    current(n) = stack % a % supercurrent(0,1)
  end do

  ! Calculate the critical current
  critical = maxval(abs(current))

  ! Write out the final results
  call finalize



  !--------------------------------------------------------------------------------!
  !                                 SUBROUTINES                                    !
  !--------------------------------------------------------------------------------!

contains
  impure subroutine prehook
    ! Write out status information.
    call status_body('Phase difference', phase(n))
  end subroutine

  impure subroutine posthook
    ! Write results to output files.
    character(len=5) :: filename
    write(filename,'(f5.3)') phase(n)
!   call stack % write_supercurrent('current.' // filename // '.dat')
!   call stack % write_gap(         'gap.'     // filename // '.dat')
  end subroutine

  impure subroutine finalize
    ! Status information
    call status_head('CRITICAL CURRENT')
    call status_body('Result', critical)
    call status_foot

    ! Write the current-phase relation to file
    call dump('current.dat', [phase, current], ['Phase             ', 'Charge current    '])

    ! Write the critical current to file
    call dump('critical.dat', critical)

    ! Close output files
    close(stdout)
    close(stderr)
  end subroutine
end program
