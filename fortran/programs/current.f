!> Author:   Jabir Ali Ouassou
!> Date:     2016-03-22
!> Category: Programs
!>
!> This program calculates the charge and spin currents in a Josephson junction as a function of the
!> phase difference between the two surrounding superconductors. In particular, the critical current
!> of the junction is also calculated. The heterostructure is constructed based on the configuration
!> file 'materials.conf', which the program expects to find in the runtime directory. The results are
!> then written to various output files that will be created in the runtime directory.

program critical_current
  use :: structure_m
  use :: stdio_m
  use :: math_m
  use :: calculus_m

  !--------------------------------------------------------------------------------!
  !                                GLOBAL VARIABLES                                !
  !--------------------------------------------------------------------------------!

  ! Declare the superconducting structure
  type(structure)                 :: stack
  type(superconductor), target    :: sa, sb

  ! Declare program control parameters
  integer,  parameter             :: iterations  = 102
  real(wp), parameter             :: tolerance   = 1e-6_wp
  real(wp), parameter             :: periodicity = 4.00_wp

  ! Declare variables used by the program
  character(len=132)              :: filename   = ''
  real(wp), dimension(iterations) :: phase      = 0
  real(wp), dimension(iterations) :: current    = 0
  real(wp)                        :: critical   = 0
  integer                         :: unit       = 0
  integer                         :: iostat     = 0
  integer                         :: n          = 0



  !--------------------------------------------------------------------------------!
  !                           INITIALIZATION PROCEDURE                             !
  !--------------------------------------------------------------------------------!

  ! Construct the central material stack
  stack = structure('materials.conf')

  ! Construct the surrounding superconductors
  sa = superconductor()
  sb = superconductor()

  ! Lock the superconductors from updates
  call sa % conf('order','-1')
  call sb % conf('order','-1')

  ! Connect the superconductors to the stack
  stack % a % material_a => sa
  stack % b % material_b => sb



  !--------------------------------------------------------------------------------!
  !                           LINEAR SEARCH PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Calculate which phase differences to check
  call linspace(phase, 1e-6_wp, periodicity-1e-6_wp)

  ! Calculate the charge current as a function of phase difference
  do n=1,iterations
    ! Status information
    call status_head('UPDATING PHASE')
    call status_body('Phase difference', phase(n))
    call status_foot

    ! Update the phase
    call sa % init( gap = exp(((0.0,-0.5)*pi)*phase(n)) )
    call sb % init( gap = exp(((0.0,+0.5)*pi)*phase(n)) )

    ! Update the state
    call stack % update
    do while (stack % difference() > tolerance)
      ! Status information
      call status_head('UPDATING STATE')
      call status_body('Phase difference', phase(n))
      call status_body('State difference', stack % difference())
      call status_foot

      ! Update the state
      call stack % update
    end do

    ! Write all currents to file
    write(filename,'(a,f5.3,a)') 'current.', phase(n), '.dat'
    call stack % write_current(filename)

    ! Write the superconducting gap to file
    write(filename,'(a,f5.3,a)') 'gap.', phase(n), '.dat'
    call stack % write_gap(filename)

    ! Save the charge current to array
    current(n) = stack % a % current(0,1)
  end do

  ! Calculate the critical current
  critical = maxval(abs(current))

  ! Status information
  call status_head('CRITICAL CURRENT')
  call status_body('Result', critical)
  call status_foot

  ! Write the current-phase relation to file
  filename = 'current.dat'
  open(newunit = unit, file = filename, iostat = iostat, action = 'write', status = 'replace')
  if (iostat /= 0) then
    call error('Failed to open output file "' // filename // '"!')
  end if
  write(unit,'(*(a,:,"	"))') '# Phase             ', '  Charge current    '
  do n=1,size(phase)
    write(unit,'(*(es20.12e3,:,"	"))') phase(n), current(n)
  end do
  close(unit = unit)

  ! Write the critical current to file
  filename = 'critical.dat'
  open(newunit = unit, file = filename, iostat = iostat, action = 'write', status = 'replace')
  if (iostat /= 0) then
    call error('Failed to open output file "' // filename // '"!')
  end if
  do n=1,command_argument_count()
    call get_command_argument(n, filename)
    write(unit,'(a,"	")',advance='no') trim(filename)
  end do
  write(unit,'(es20.12e3)') critical
  close(unit = unit)
end program
