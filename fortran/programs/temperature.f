! This program calculates the critical temperature in a multilayer structure.
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2016-03-09
! Updated: 2016-03-23

program critical_temperature
  use mod_structure
  use mod_math
  implicit none

  !--------------------------------------------------------------------------------!
  !                                GLOBAL VARIABLES                                !
  !--------------------------------------------------------------------------------!

  ! Declare the superconducting structure
  type(structure) :: stack

  ! Declare program control parameters
  integer         :: information  = 0
  integer         :: bisections   = 12
  integer         :: iterations   = 6
  real(wp)        :: initgap      = 1e-5_wp
  real(wp)        :: minimum      = 0.00_wp
  real(wp)        :: maximum      = 1.00_wp

  ! Declare variables used by the program
  real(wp)        :: temperature
  integer         :: n



  !--------------------------------------------------------------------------------!
  !                          INITIALIZATION PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Construct the multilayer stack based on a config file
  stack = structure('structure.conf')

  ! Initialize the stack to a barely superconducting state
  call stack % init(cx(initgap,0.0_wp))

  ! Bootstrap the materials at zero temperature
  do while(stack % difference() > 1e-7)
    call stack % update(freeze = .true.)
  end do

  ! Save the current state of the materials
  call stack % save

  !--------------------------------------------------------------------------------!
  !                           BINARY SEARCH PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  do n = 1,bisections
    ! Load the saved material states
    call stack % load

    ! Set the temperature
    temperature = (minimum+maximum)/2
    call stack % temperature(temperature)
    write(*,*) temperature

    ! Update the stack
    call stack % update

    minimum = temperature
  end do

  ! Set the system temperature to the midpoint of the search space
  !s(1) % temperature = (minimum + maximum)/2.0_wp

  !! Perform the binary search
  !do n = 1,bisections
  !  ! Status information
  !  call print_status('           CONVERGENCE', bisection=n, iteration=0, temperature = s(1)%temperature)

  !  ! Load the states from backup
  !  if (size(f) > 0) then
  !    call s(1) % load
  !    do m = 1,size(f)
  !      call f(m) % load
  !    end do
  !  else
  !    call s(1) % init( gap = cmplx(initgap,0,kind=wp) )
  !  end if

  !  ! Update the superconductor once
  !  call s(1) % update

  !  ! Perform any additional iterations over the structure
  !  if (size(f) > 0) then
  !    do m = 1,iterations
  !      ! Status information
  !      call print_status('           CONVERGENCE', bisection=n, iteration=m, temperature = s(1)%temperature)

  !      ! Update the ferromagnets right-to-left (including edges)
  !      do k=size(f),1,-1
  !        call f(k) % update
  !      end do

  !      ! Update the ferromagnets left-to-right (excluding edges)
  !      if (size(f) > 1) then
  !        do k=2,size(f)-1
  !          call f(k) % update
  !        end do
  !      end if

  !      ! Update the superconductor
  !      call s(1) % update
  !    end do
  !  end if

  !  ! Check whether the gap has increased or not, and update the search space accordingly
  !  if (abs(s(1)%get_gap_mean()/initgap) >= 1.0_wp) then
  !    minimum = s(1) % temperature
  !  else
  !    maximum = s(1) % temperature
  !  end if

  !  ! Set the system temperature to the miwpoint of the updated search space
  !  s(1) % temperature = (maximum + minimum)/2.0_wp
  !end do

  !! Print the final results
  !call print_status('            CONVERGED', temperature=s(1)%temperature)
end program
