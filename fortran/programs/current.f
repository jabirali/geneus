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
  real(wp), parameter             :: threshold   = 1e-2
  real(wp), parameter             :: tolerance   = 1e-6
  integer,  parameter             :: iterations  = 51

  ! Declare variables used by the program
  character(len=132)                  :: filename = ''
  real(wp), dimension(:), allocatable :: phase
  real(wp), dimension(:), allocatable :: current
  real(wp)                            :: critical
  integer                             :: materials
  integer                             :: superconductors
  logical                             :: bootstrap
  logical                             :: shortcircuit
  integer                             :: n



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

  ! Count the number of materials
  materials       = stack % materials()
  superconductors = stack % superconductors()

  ! If we try to solve the diffusion and selfconsistency equations simultaneously,
  ! using bad initial guesses for the order parameters, the numerics can become 
  ! unstable and converge slower. This is especially problematic in systems with
  ! multiple superconductors, where this may result in the numerics diverging.
  ! It can therefore be useful to start with a bootstrap procedure, where we
  ! solve the diffusion equation for fixed superconducting gaps, before doing
  ! a fully selfconsistent calculation. Here, we check if that is necessary.
  bootstrap = (superconductors > 0)

  ! If we have an unlocked superconductor in the system, then multiple iterations
  ! are required to find a simultaneous solution of the diffusion equation and the
  ! selfconsistency equation. If we have at least two materials in the system, then
  ! multiple iterations are required to make the diffusion equation solutions in
  ! each material consistent with each other. Here, we check if that's unnecessary.
  shortcircuit = (materials == 1) .and. (superconductors == 0)

  ! Depending on the number of superconductors in the junction, we might get a
  ! 2πn-periodicity instead of a 2π-periodicity, and should account for that.
  allocate(phase((superconductors+1)*(iterations-1)+1))
  allocate(current(size(phase)))
  call linspace(phase, 1e-6_wp, 2*(superconductors+1) - 1e-6_wp)



  !--------------------------------------------------------------------------------!
  !                              BOOTSTRAP PROCEDURE                               !
  !--------------------------------------------------------------------------------!

  ! Solve the diffusion equations
  if (bootstrap) then
    do
      ! Status information
      call status_head('BOOTSTRAPPING')
      call status_body('State difference', stack % difference())
      call status_foot

      ! Update the material state (non-selfconsistently)
      call stack % update(freeze = .true.)

      ! Stop if all materials have converged within the threshold
      if (stack % difference() < threshold) then
        exit
      end if
    end do
  end if

  ! Solve the selfconsistency equations (if necessary)
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
    do
      ! Status information
      call status_head('UPDATING STATE')
      call status_body('Phase difference', phase(n))
      call status_body('State difference', stack % difference())
      call status_foot

      ! Update the state
      call stack % update

      ! Stop if all materials have converged within the tolerance
      if (shortcircuit .or. stack % difference() < tolerance) then
        exit
      end if
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
