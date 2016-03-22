!! This program calculates the critical current in an S/X/S heterostructure as a function of the phase
!! difference between the two superconductors. The heterostructure is constructed based in the config
!! file 'simulation.conf', and the output is written to the files 'current.dat' and 'critical.dat'.
!!
!! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
!! Created: 2016-03-22
!! Updated: 2016-03-22

program critical_current
  use mod_structure
  use mod_stdio
  use mod_math
  implicit none

  !--------------------------------------------------------------------------------!
  !                                GLOBAL VARIABLES                                !
  !--------------------------------------------------------------------------------!

  ! Declare the superconducting structure
  type(structure)                     :: stack
  type(superconductor), target        :: sa, sb

  ! Declare program control parameters
  integer,  parameter                 :: iterations = 101
  real(wp), parameter                 :: tolerance  = 1d-6

  ! Declare variables used by the program
  real(wp), dimension(iterations)     :: current     = 0
  real(wp)                            :: critical    = 0
  real(wp), dimension(100*iterations) :: domain
  real(wp), dimension(iterations)     :: phase
  integer                             :: n



  !--------------------------------------------------------------------------------!
  !                          INITIALIZATION PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Construct the multilayer stack based on a config file
  stack = structure('simulation.conf')

  ! Construct two superconductors
  sa = superconductor(30.0_wp)
  sb = superconductor(30.0_wp)

  ! Lock the new superconductors
  sa % lock = .true.
  sb % lock = .true.

  ! Connect the superconductors to the stack
  stack % a % material_a => sa
  stack % b % material_b => sb



  !--------------------------------------------------------------------------------!
  !                           LINEAR SEARCH PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Calculate which phase differences to check
  call linspace(phase,  1d-6, 1-1d-6)
  call linspace(domain, 1d-6, 1-1d-6)

  ! Calculate the charge current as a function of phase difference
  do n=1,iterations
    ! Update the phase
    call sa % init( gap = exp(((0.0,-0.5)*pi)*phase(n)) )
    call sb % init( gap = exp(((0.0,+0.5)*pi)*phase(n)) )

    ! Update the system
    call stack % update
    do while (stack % difference() > tolerance)
      call stack % update
    end do

    ! Calculate the current
    current(n) = stack % a % current(0,1)
    write(*,*) phase(n), current(n)

    ! Flush output
    flush(stdout)
  end do

  ! Calculate the critical current
  critical = maxval(abs(current))
  write(*,*) 'Critical current:', critical

  ! Interpolate the critical current
  critical = maxval(abs(interpolate(phase, current, domain)))
  write(*,*) 'Interpolated:', critical

  ! Write current(0:3,:) to current.dat, and critical to critical.dat
end program
