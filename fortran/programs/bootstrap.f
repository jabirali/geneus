!> Author:   Jabir Ali Ouassou
!> Date:     2016-03-08
!> Category: Programs
!>
!> This program calculates various physical observables for a superconducting thin-film structure
!> in equilibrium, including the charge and spin currents, density of states, and order parameter.
!> The heterostructure is constructed based on the configuration file 'materials.conf',  which is
!> expected to be in the runtime directory. The output is written to files in the same directory.
!> In contrast to equilibrium.f, this version of the program does a non-selfconsistent calculation
!> first, and then continues with a fully selfconsistent solution of the problem after convergence.

program bootstrap
  use :: structure_m
  use :: stdio_m

  ! Create the superconducting structure
  type(structure) :: stack
  stack = structure('materials.conf')

  ! Non-selfconsistent solution
  do
    ! Status information
    call status_head('BOOTSTRAPPING')
    call status_body('State difference', stack % difference())
    call status_foot

    ! Update the material state
    call stack % update(freeze = .true.)

    ! Stop if we have weak convergence
    if (stack % difference() < 1e-2) then
      exit
    end if
  end do

  ! Selfconsistent solution
  do
    ! Status information
    call status_head('UPDATING STATE')
    call status_body('State difference', stack % difference())
    call status_foot

    ! Update the material state
    call stack % update

    ! Write the results to files
    call stack % write_density('density.dat')
    call stack % write_current('current.dat')
    call stack % write_magnetization('magnetization.dat')
    call stack % write_gap('gap.dat')

    ! Stop if we have convergence
    if (stack % difference() < 1e-9) then
      exit
    end if
  end do

  ! Status information
  call status_head('EQUILIBRIUM')
  call status_body('State difference', stack % difference())
  call status_foot
end program
