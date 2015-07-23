! CR: 2015-07-23
! UP   -- " --

module mod_critical
  use mod_system
  use mod_superconductor
  use mod_multilayer
  implicit none
contains

subroutine critical_temperature(material, bisections, iterations, lower, upper, gap)
  ! This subroutine calculates the critical temperature of a superconducting material using a binary search algorithm.
  ! The two compulsory arguments are 'material', which is the superconducting material that we will find the critical
  ! temperature of, and 'update', which is a subroutine that updates the current system state by solving the equations
  ! of motion in the entire hybrid system. For more fine-grained control over the process, the user can also specify
  ! the number of temperatures to check during the binary search (bisections), the number of updates to perform per
  ! temperature until the state stabilizes (iterations), and what temperature range to investigate (lower to upper).
  ! Finally, 'gap' determines the superconducting gap that is used to reinitialize the materials at all temperatures.

  ! Declare data types for input variables
  class(superconductor), target,   intent(inout) :: material
  integer,               optional, intent(in   ) :: bisections
  integer,               optional, intent(in   ) :: iterations
  real(dp),              optional, intent(in   ) :: lower 
  real(dp),              optional, intent(in   ) :: upper
  complex(dp),           optional, intent(in   ) :: gap

  ! Declare data types for internal variables
  integer                                        :: bisections_ = 12
  integer                                        :: iterations_ =  5
  real(dp)                                       :: lower_      =  0.0000_dp
  real(dp)                                       :: upper_      =  1.0000_dp
  complex(dp)                                    :: gap_        =  0.0001_dp
  integer                                        :: n, m

  ! Handle optional arguments
  call process_arguments

  ! Set the system temperature to the midpoint of the search space
  call material%set_temperature( (upper_ + lower_)/2.0_dp )

  ! Perform the binary search
  do n = 1,bisections_
    ! Status information
    call print_information

    ! Initialize a weakly superconducting state
    call initialize_all(material, gap_)

    ! Update the state of the multilayer structure
    do m = 1,iterations_
      call update_all(material)
    end do

    ! Check whether the gap has increased or not, and update the search space accordingly
    if (abs(material%get_gap_mean()/gap_) >= 1.0_dp) then
      lower_ = material%get_temperature()
    else
      upper_ = material%get_temperature()
    end if

    ! Set the system temperature to the midpoint of the updated search space
    call material%set_temperature( (upper_ + lower_)/2.0_dp )
  end do

  ! Print the final results
  call print_results
contains
  subroutine process_arguments
    ! This nested subroutine processes the optional arguments to the parent subroutine.
    if (present(bisections)) then
      bisections_ = bisections
    end if

    if (present(iterations)) then
      iterations_ = iterations
    end if

    if (present(lower)) then
      lower_ = lower
    end if

    if (present(upper)) then
      upper_ = upper
    end if

    if (present(gap)) then
      gap_ = gap
    end if
  end subroutine

  subroutine print_information
    ! Determine how much CPU time has elapsed
    real(sp) :: time
    call cpu_time(time)

    ! Print the progress information to standard out
    write(*,'(a)') '                                     '
    write(*,'(a)') '╒═══════════════════════════════════╕'
    write(*,'(a)') '│       PROGRESS  INFORMATION       │'
    write(*,'(a)') '├───────────────────────────────────┤'
    write(*,'(a,6x,a,i2.2,a,i2.2,7x,a)')                &
      '│','Binary search: ', n, ' / ', bisections,     '│'
    write(*,'(a,6x,a,f8.6,6x,a)')                       &
      '│','Temperature:   ',material%get_temperature(),'│'
    write(*,'(a,6x,a,i2.2,a,i2.2,a,i2.2,6x,a)')         &
      '│','Elapsed time:  ',                            &
      int(time/3600.0_sp),':',                          &
      int(mod(time,3600.0_sp)/60.0_sp),':',             &
      int(mod(time,60.0_sp)),                          '│'
    write(*,'(a)') '╘═══════════════════════════════════╛'
  end subroutine

  subroutine print_results
    write(*,'(a)') '                                     '
    write(*,'(a)') '╒═══════════════════════════════════╕'
    write(*,'(a)') '│           FINAL RESULTS           │'
    write(*,'(a)') '├───────────────────────────────────┤'
    write(*,'(a)') '│       Critical temperature:       │'
    write(*,'(a,8x,f18.16,9x,a)') '│', material%get_temperature(), '│'
    write(*,'(a)') '╘═══════════════════════════════════╛'
  end subroutine
end subroutine

end module
