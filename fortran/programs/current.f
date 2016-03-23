!! This program calculates the charge and spin currents in an Josephson junction as a function of the
!! phase difference between the two surrounding superconductors.  In particular, the critical current
!! of the junction is calculated.  The heterostructure is constructed based on the configuration file
!! 'structure.conf', which the program expects to find in the runtime directory. The results are then
!! written to various output files that will be created in the same directory.
!!
!! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
!! Created: 2016-03-22
!! Updated: 2016-03-23

program critical_current
  use mod_structure
  use mod_stdio
  use mod_math
  implicit none

  !--------------------------------------------------------------------------------!
  !                                GLOBAL VARIABLES                                !
  !--------------------------------------------------------------------------------!

  ! Declare the superconducting structure
  type(structure)                 :: stack
  type(superconductor), target    :: sa, sb

  ! Declare program control parameters
  integer,  parameter             :: iterations = 21
  real(wp), parameter             :: tolerance  = 1e-6
  logical                         :: loop       = .true.

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
  stack = structure('structure.conf')

  ! Construct the surrounding superconductors
  sa = superconductor()
  sb = superconductor()

  ! Lock the superconductors from updates
  call sa % conf('lock','T')
  call sb % conf('lock','T')

  ! Connect the superconductors to the stack
  stack % a % material_a => sa
  stack % b % material_b => sb



  !--------------------------------------------------------------------------------!
  !                           LINEAR SEARCH PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Calculate which phase differences to check
  call linspace(phase, 1d-6, 1-1d-6)

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
    if (loop) then
      do while (stack % difference() > tolerance)
        ! Status information
        call status_head('UPDATING STATE')
        call status_body('Phase difference', phase(n))
        call status_body('State difference', stack % difference())
        call status_foot

        ! Update the state
        call stack % update
      end do
    end if

    ! Write all currents to file
    write(filename,'(a,f5.3,a)') 'current.', phase(n), '.dat'
    call stack % write_current(filename)

    ! Save the charge current
    current(n) = stack % a % current(0,1)
  end do

  ! Calculate the critical current
  critical = maxval(abs(current))

  ! Status information
  call status_head('COMPLETE')
  call status_body('Critical current', critical)
  call status_foot

  ! Write the current-phase relation to file
  filename = 'current.dat'
  open(newunit = unit, file = filename, iostat = iostat, action = 'write', status = 'replace')
  if (iostat /= 0) then
    call error('Failed to open output file "' // filename // '"!')
  end if
  do n=1,size(phase)
    write(unit,*) phase(n), current(n)
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
    write(unit,'(a,2x)',advance='no') trim(filename)
  end do
  write(unit,*) critical
  close(unit = unit)
end program
