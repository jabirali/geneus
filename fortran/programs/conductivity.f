!> Author:   Jabir Ali Ouassou
!> Date:     2017-02-15
!> Category: Programs
!>
!> This program calculates the charge conductivity of an S/X/N superconducting thin-film structure.

program conductivity
  use :: structure_m
  use :: stdio_m
  use :: math_m

  ! Declare the superconducting structure
  type(structure) :: stack

  ! Declare program control parameters
  real(wp), parameter :: threshold = 1e-2
  real(wp), parameter :: tolerance = 1e-8

  ! Redefine stdout and stderr 
  stdout = output('output.log')
  stderr = output('error.log')

  ! Construct the superconducting structure
  stack = structure('materials.conf')

  ! @TODO: Check that the stack has only one layer

  ! @TODO: Check if the boundary conditions are spin-active?

  ! @TODO: Connect S and N layers around it

  ! @TODO: Disable updates in the S and N

  ! @TODO: Make the N interface transparent

  ! Non-selfconsistent bootstrap procedure
  call stack % converge(threshold = threshold, bootstrap = .true.)

  ! Selfconsistent convergence procedure
  call stack % converge(threshold = tolerance, posthook = posthook)

  ! Write out the final results
  call finalize

  ! @TODO: Calculate the conductivity
  !
  !    Trapezoid integration:
  !     integral[f(z), z=a, z=b] = 0.5·sum( (z(k+1)-z(k))·(f(k+1)-f(k)) )
  !    
  !    Implementation:
  !     h0(0)  = h_alpha
  !     h0(1:) = h(0) + 0.5*cumsum( (z(1:n)-z(0:n-1)) * (Minv(1:n)-Minv(0:n-1)) )
  !    
  !     h(0)  = h0(0)
  !     h(1:) = h0(1:) + i[A, 0.5*cumsum( (z(1:n)-z(0:n-1)) * (h(1:n)-h(0:n-1)) )]



  !--------------------------------------------------------------------------------!
  !                                 SUBROUTINES                                    !
  !--------------------------------------------------------------------------------!

contains
  impure subroutine posthook
    ! Write results to output files.
    use :: structure_m

    call stack % write_density('density.dat')
  end subroutine

  impure subroutine finalize
    ! Write out the final results.
    use :: stdio_m

    ! Status information
    call status_head('COMPLETE')
    call status_body('State difference', stack % difference())
    call status_body('Charge violation', stack % chargeviolation())
    call status_foot

    ! Close output files
    close(stdout)
    close(stderr)
  end subroutine
end program
