program equilibrium
  use :: structure_m
  use :: stdio_m

  ! Declare variables
  integer        :: n
  character(128) :: tmp

  ! Create the superconducting structure
  type(structure) :: stack

  call stack % push('superconductor')
  call stack % conf('length', '1.00')
  call stack % conf('conductance_b', '0.40')
  call stack % conf('polarization_b', '0.999')
  call stack % conf('magnetization_b', '[0,0,1]')
  call stack % conf('misalignment_b', '[sin(0.5*pi),0,cos(0.5*pi)]')

  call stack % push('halfmetal')
  call stack % conf('length', '12.00')
  call stack % conf('polarization', '0.999')
  call stack % conf('conductance_a', '0.40')

  do n=0,10000
    ! Update the interface parameters
    write(tmp,*) 2.50 * (n/10000.0) + 1e-12
    call stack % a % conf('spinmixing_b', tmp)
    call stack % b % conf('spinmixing_a', tmp)

    ! Status information
    call status_head('PARAMETER UPDATE')
    call status_body('Spinmixing', 2.50 * (n/10000.0))
    call status_foot

    ! Main loop
    do
      ! Status information
      call status_head('STATE UPDATE')
      call status_body('State difference', stack % difference())
      call status_foot

      ! Update the material state
      call stack % update

      ! Write the results to files
      call stack % write_density('density.dat')
      call stack % write_current('current.dat')
      call stack % write_gap('gap.dat')

      ! Stop if we have convergence
      if (stack % difference() < 1e-3) then
        exit
      end if
    end do
  end do

  ! Strong convergence
  do
    ! Status information
    call status_head('FINAL UPDATE')
    call status_body('State difference', stack % difference())
    call status_foot

    ! Update the material state
    call stack % update

    ! Write the results to files
    call stack % write_density('density.dat')
    call stack % write_current('current.dat')
    call stack % write_gap('gap.dat')

    ! Stop if we have convergence
    if (stack % difference() < 1e-6) then
      exit
    end if
  end do

  ! Status information
  call status_head('CONVERGED')
  call status_body('State difference', stack % difference())
  call status_foot
end program
