!> Author:   Jabir Ali Ouassou
!> Date:     2016-03-24
!> Category: Programs
!>
!> This program calculates the critical temperature of an arbitrary superconducting hybrid structure
!> for a given set of physical parameters. In order to obtain the critical temperature as a function
!> of these parameters, the program has to be invoked multiple times with different input parameters.
!> The structure is constructed based on the configuration file 'materials.conf', which the program
!> expects to find in the runtime directory, and the results are written to the file 'critical.dat'. 

program critical_temperature
  use structure_m
  use stdio_m
  use math_m
  implicit none

  !--------------------------------------------------------------------------------!
  !                                GLOBAL VARIABLES                                !
  !--------------------------------------------------------------------------------!

  ! Declare the superconducting structure
  type(structure)                 :: stack

  ! Declare program control parameters
  character(*), parameter         :: filename   = 'critical.dat'
  integer,      parameter         :: bootstraps = 12
  integer,      parameter         :: bisections = 12
  integer,      parameter         :: iterations = 06
  real(wp),     parameter         :: tolerance  = 1e-8_wp
  real(wp),     parameter         :: initgap    = 1e-5_wp

  ! Declare variables used by the program
  character(len=132)              :: argument   = ''
  real(wp)                        :: minimum    = 0.00_wp
  real(wp)                        :: maximum    = 1.00_wp
  real(wp)                        :: critical   = 0.50_wp
  integer                         :: unit       = 0
  integer                         :: iostat     = 0
  integer                         :: n          = 0
  integer                         :: m          = 0



  !--------------------------------------------------------------------------------!
  !                           INITIALIZATION PROCEDURE                             !
  !--------------------------------------------------------------------------------!

  ! Construct the material stack
  stack = structure('materials.conf')

  ! Initialize the stack to a barely superconducting state
  call stack % init(cx(initgap,0.0_wp))

  ! Check the number of materials in the stack
  m = stack % count()

  ! Bootstrap the material states at zero temperature
  do
    ! Update counter
    n = n+1

    ! Status information
    call status_head('INITIALIZING')
    call status_body('Temperature', 0.0)
    call status_body('State difference', stack % difference())
    call status_body('Iteration', n)
    call status_foot

    ! Update materials
    call stack % update(freeze = .true.)

    ! If we only have one layer, then one update is sufficient
    if (m <= 1) then
      exit
    end if

    ! If we have multiple layers, check the iteration number
    if (n >= bootstraps) then
      ! Check for convergence
      if (stack % difference() < tolerance) then
        exit
      end if
    end if
  end do

  ! Save the current state of the materials
  call stack % save



  !--------------------------------------------------------------------------------!
  !                           BINARY SEARCH PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  do n = 1,bisections
    ! Set the temperature of the materials
    call stack % temperature(critical)

    ! Load the saved material states
    call stack % load

    ! Update the material states
    do m = 1,iterations
      ! Status information
      call status_head('UPDATING STATE')
      call status_body('Temperature', critical)
      call status_body('Bisection',   n)
      call status_body('Iteration',   m)
      call status_foot

      ! Update the stack
      call stack % update
    end do

    ! Update the critical temperature estimate
    if (stack % gap() >= initgap) then
      minimum = critical
    else
      maximum = critical
    end if
    critical = (minimum + maximum)/2
  end do

  ! Status information
  call status_head('CRITICAL TEMPERATURE')
  call status_body('Result', critical)
  call status_foot

  ! Write the critical temperature to file
  open(newunit = unit, file = filename, iostat = iostat, action = 'write', status = 'replace')
  if (iostat /= 0) then
    call error('Failed to open output file "' // filename // '"!')
  end if
  do n=1,command_argument_count()
    call get_command_argument(n, argument)
    write(unit,'(a,"	")',advance='no') trim(argument)
  end do
  write(unit,'(es20.12e3)') critical
  close(unit = unit)
end program
