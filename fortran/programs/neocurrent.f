!! This program calculates the critical current in an S/X/S heterostructure as a function of the phase
!! difference between the two superconductors. The heterostructure is constructed based in the config
!! file 'simulation.conf', and the output is written to the file 'critical_current.dat'.
!!
!! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
!! Created: 2016-03-22
!! Updated: 2016-03-22

program critical_current
  use mod_structure
  use mod_superconductor
  use mod_stdio
  use mod_math

  ! Declare variables 
  type(structure) :: stack
  real(wp)        :: phase
  integer         :: n

  ! Construct the multilayer stack based on a config file
  stack = structure('simulation.conf')

  ! Verify that this is an S/X/S junction of some kind, and lock the superconductors
  associate(a => stack % a, b => stack % b)
    if (associated(stack % a, stack % b)) then
      call error('Minimum two superconductors required for the calculations!')
    end if

    select type(a)
      class is (superconductor)
        call a % conf('lock', 'T')
      class default
        call error('First material in the heterostructure is not a superconductor!')
    end select

    select type(b)
      class is (superconductor)
        call b % conf('lock', 'T')
      class default
        call error('Last material in the heterostructure is not a superconductor!')
    end select

    if (associated(a % material_b, b)) then
      call error('Minimum one material required between the two superconductors!')
    end if
  end associate

  do n=0,100
    ! Update the phase
    phase = (n+1e-6)/(100+2e-6)
    call stack % a % init( gap = exp(((0.0,-0.5)*pi)*phase) )
    call stack % b % init( gap = exp(((0.0,+0.5)*pi)*phase) )

    ! Update the system
    call stack % update
    do while (stack % difference() > 1e-4)
      call stack % update
    end do

    ! Calculate the critical current
    write(*,*) phase, abs(stack % a % material_b % current(1,1))

    ! Flush output
    flush(stdout)
  end do
end program
