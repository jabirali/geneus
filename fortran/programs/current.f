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
  integer,  parameter             :: iterations  = 51
  real(wp), parameter             :: tolerance   = 1e-6_wp

  ! Declare variables used by the program
  character(len=132)                  :: filename = ''
  real(wp), dimension(:), allocatable :: phase
  real(wp), dimension(:), allocatable :: current
  real(wp)                            :: critical
  integer                             :: n
  integer                             :: m



  !--------------------------------------------------------------------------------!
  !                           INITIALIZATION PROCEDURE                             !
  !--------------------------------------------------------------------------------!

  ! Construct the central material stack
  stack = structure('materials.conf')

  ! Construct the surrounding superconductors
  sa = superconductor()
  sb = superconductor()

  ! Lock the superconductors from updates
  call sa % conf('order','0')
  call sb % conf('order','0')

  ! Connect the superconductors to the stack
  stack % a % material_a => sa
  stack % b % material_b => sb

  ! Check the number of unlocked superconductors. This is used to determine the
  ! maximum periodicity the current-phase relation can have, and therefore which
  ! phase-differences we should check when solving the diffusion equations.
  m = stack % superconductors()
  allocate(phase((m+1)*(iterations-1)+1))
  allocate(current(size(phase)))
  call linspace(phase, 1e-6_wp, 2*(m+1) - 1e-6_wp)

  ! Solve the selfconsistency equations
  call stack % load()


  !--------------------------------------------------------------------------------!
  !                           LINEAR SEARCH PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Calculate the charge current as a function of phase difference
  do n=1,size(phase)
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



  !--------------------------------------------------------------------------------!
  !                            FINALIZATION PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Calculate the critical current
  critical = maxval(abs(current))

  ! Status information
  call status_head('CRITICAL CURRENT')
  call status_body('Result', critical)
  call status_foot

  ! Write the current-phase relation to file
  call dump('current.dat', [phase, current], ['Phase             ', 'Charge current    '])

  ! Write the critical current to file
  call dump('critical.dat', critical)

  ! Deallocate memory
  deallocate(phase, current)
end program
