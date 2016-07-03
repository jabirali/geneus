program equilibrium
  use :: structure_m
  use :: stdio_m

  ! Declare variables
  integer        :: n
  character(128) :: tmp

  ! Create the superconducting structure
  type(structure) :: stack

  call stack % push('superconductor')
  call stack % conf('length', '3.00')
  call stack % conf('spinmixing_b', '1.00')
  call stack % conf('polarization_b', '0.999')
  call stack % conf('magnetization_b', '[0,0,1]')
  call stack % conf('misalignment_b', '[sin(0.5*pi),0,cos(0.5*pi)]')

  call stack % push('halfmetal')
  !call stack % conf('length', '20.00')
  !call stack % conf('polarization', '0.999')

  call stack % update

  !do n=1,2
  !  !write(tmp,*) n/1000.0
  !  write(*,*) tmp
  !  call stack % a % conf('conductance_b', tmp)
  !end do

  ! Main loop
  !do
  !  ! Status information
  !  call status_head('SELFCONSISTENT UPDATE')
  !  call status_body('State difference', stack % difference())
  !  call status_foot

  !  ! Update the material state
  !  call stack % update

  !  ! Write the results to files
  !  call stack % write_density('density.dat')
  !  call stack % write_current('current.dat')
  !  call stack % write_gap('gap.dat')

  !  ! Stop if we have convergence
  !  if (stack % difference() < 1e-6) then
  !    exit
  !  end if
  !end do

  !! Status information
  !call status_head('CONVERGED')
  !call status_body('State difference', stack % difference())
  !call status_foot
end program
