!! This program calculates the charge and spin currents in an S/X/S structure as a function of the
!! phase difference between the two superconductors, and in particular the critical current of the
!! junction. This heterostructure is constructed based on the configuration file 'structure.conf',
!! and the output is written to a series of output files in the same directory.
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
  integer,  parameter                 :: iterations = 6
  real(wp), parameter                 :: tolerance  = 1d-6

  ! Declare variables used by the program
  character(len=132)                  :: filename    = ''
  real(wp)                            :: critical    = 0
  real(wp), dimension(iterations)     :: current     = 0
  real(wp), dimension(100*iterations) :: domain
  real(wp), dimension(iterations)     :: phase
  integer                             :: n



  !--------------------------------------------------------------------------------!
  !                          INITIALIZATION PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Construct the central material stack
  stack = structure('structure.conf')

  ! Construct the surrounding superconductors
  sa = superconductor()
  sb = superconductor()

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

    ! Write results to file
    write(filename,'(a,f5.3,a)') 'current.', phase(n), '.dat'
    call stack % write_current(filename)

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
