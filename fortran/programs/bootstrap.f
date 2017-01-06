!> Author:   Jabir Ali Ouassou
!> Date:     2016-03-08
!> Category: Programs
!>
!> This program calculates various physical observables for a superconducting thin-film structure
!> in equilibrium, including the charge and spin currents, density of states, and order parameter.
!> The heterostructure is constructed based on the configuration file 'materials.conf',  which is
!> expected to be in the runtime directory. The output is written to files in the same directory.

program equilibrium
  use :: structure_m
  use :: stdio_m
  use :: math_m

  !--------------------------------------------------------------------------------!
  !                                 INITIALIZATION                                 !
  !--------------------------------------------------------------------------------!

  ! Declare program control parameters
  real(wp), parameter :: threshold = 1e-2
  real(wp), parameter :: tolerance = 1e-6

  ! Declare variables used by the program
  type(structure)     :: stack
  integer             :: materials
  integer             :: superconductors
  logical             :: bootstrap
  logical             :: shortcircuit

  ! Create the superconducting structure
  stack = structure('materials.conf')

  ! Count the number of materials
  materials       = stack % materials()
  superconductors = stack % superconductors()

  ! If we try to solve the diffusion and selfconsistency equations simultaneously,
  ! using bad initial guesses for the order parameters, the numerics can become 
  ! unstable and converge slower. This is especially problematic in systems with
  ! multiple superconductors, where this may result in the numerics diverging.
  ! It can therefore be useful to start with a bootstrap procedure, where we
  ! solve the diffusion equation for fixed superconducting gaps, before doing
  ! a fully selfconsistent calculation. This is only necessary in systems with
  ! at least two materials, where at least one of them is a superconductor.
  bootstrap = (materials > 1) .and. (superconductors > 0)

  ! If we have an unlocked superconductor in the system, then multiple iterations
  ! are required to find a simultaneous solution of the diffusion equation and the
  ! selfconsistency equation. If we have at least two materials in the system, then
  ! multiple iterations are required to make the diffusion equation solutions in
  ! each material consistent with each other. Here, we check if that's unnecessary.
  shortcircuit = (materials == 1) .and. (superconductors == 0)



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

  ! Solve the selfconsistency equations
  call stack % load()



  !--------------------------------------------------------------------------------!
  !                              SOLUTION PROCEDURE                                !
  !--------------------------------------------------------------------------------!

  ! Solve the diffusion and selfconsistency equations
  do
    ! Status information
    call status_head('UPDATING STATE')
    call status_body('State difference', stack % difference())
    call status_foot

    ! Update the material state (selfconsistently)
    call stack % update

    ! Write the results to files
    call stack % write_density('density.dat')
    call stack % write_current('current.dat')
    call stack % write_magnetization('magnetization.dat')
    call stack % write_gap('gap.dat')

    ! Stop if all materials have converged within the tolerance
    if (shortcircuit .or. stack % difference() < tolerance) then
      exit
    end if
  end do



  !--------------------------------------------------------------------------------!
  !                            FINALIZATION PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Status information
  call status_head('EQUILIBRIUM')
  call status_body('State difference', stack % difference())
  call status_foot
end program
