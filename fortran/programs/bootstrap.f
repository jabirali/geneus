!> Author:   Jabir Ali Ouassou
!> Date:     2016-04-14
!> Category: Programs
!>
!> This program is a version of `equilibrium.f` that performs a two-step solution of the diffusion
!> equations.  First, it solves the non-selfconsistent equations until course-grained convergence,
!> and then it proceeds to solve the full selfconsistent equations until fine-grained convergence.
!> This makes the program slower but also more robust compared to its one-step counterpart solver.

program bootstrap
  use :: structure_m
  use :: stdio_m

  ! Create the superconducting structure
  type(structure) :: stack
  stack = structure('materials.conf')

  ! Solve the non-selfconsistent problem
  do
    ! Status information
    call status_head('NON-SELFCONSISTENT UPDATE')
    call status_body('State difference', stack % difference())
    call status_foot

    ! Update materials
    call stack % update(freeze = .true.)

    ! Check if we have convergence
    if (stack % difference() < 1e-3) then
      ! Easiest way to call all update posthooks
      call stack % save
      call stack % load

      ! Write the results to files
      call stack % write_density('density.dat')
      call stack % write_current('current.dat')
      call stack % write_gap('gap.dat')

      ! Continue with the selfconsistent problem
      exit
    end if
  end do

  ! Solve the full selfconsistent problem
  do
    ! Status information
    call status_head('SELFCONSISTENT UPDATE')
    call status_body('State difference', stack % difference())
    call status_foot

    ! Update the material state
    call stack % update

    ! Write the results to files
    call stack % write_density('density.dat')
    call stack % write_current('current.dat')
    call stack % write_gap('gap.dat')

    ! Check if we have convergence
    if (stack % difference() < 1e-6) then
      exit
    end if
  end do

  ! Status information
  call status_head('CONVERGED')
  call status_body('State difference', stack % difference())
  call status_foot
end program
