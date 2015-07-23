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
    class(superconductor), target, intent(inout) :: material
    integer,     optional, intent(in   ) :: bisections
    integer,     optional, intent(in   ) :: iterations
    real(dp),    optional, intent(in   ) :: lower 
    real(dp),    optional, intent(in   ) :: upper
    complex(dp), optional, intent(in   ) :: gap

    ! Declare interfaces for input subroutines
    interface
      subroutine update
      end subroutine
    end interface

    ! Declare data types for internal variables
    integer     :: bisections_
    integer     :: iterations_
    real(dp)    :: lower_
    real(dp)    :: upper_
    complex(dp) :: gap_
    integer     :: n, m

    ! Handle optional arguments
    call process_arguments

    do n = 1,iterations_
      ! Initialize a weakly superconducting state
    end do



    !! Perform the binary search for the critical temperature
    !do n=1,iterations_
    !  ! Initialize a weakly superconducting state
    !  call material%initialize( gap_ )
    !  call material%set_gap( gap_ )

    !  ! Print status information
    !  call print_information

    !  ! Update the state of the superconductor
    !  do m=1,stabilization_
    !    write(*,'(a,i0,a,i0,a)') ' :: Updating superconductor [',m,'/',stabilization, ']'
         !call update(material)
    !    call s%update
    !  end do

    !  ! Check whether the mean gap has increased, and update the temperature bounds accordingly
    !  if (abs(s%get_gap_mean()/gap) >= 1.0_dp) then
    !    lower = s%get_temperature()
    !  else
    !    upper = s%get_temperature()
    !  end if

    !  ! Update the superconductor temperature based on the new bounds
    !  call s%set_temperature( (upper_+lower_)/2.0_dp )
    !end do

    !! Print final results
    !call print_results

    call update
    
  contains
    subroutine process_arguments
      ! .....
      if (present(bisections)) then
        bisections_ = bisections
      else
        bisections_ = 12
      end if

      if (present(iterations)) then
        iterations_ = iterations
      else
        iterations_ = 5
      end if

      if (present(lower)) then
        lower_ = lower
      else
        lower_ = 0.0_dp
      end if

      if (present(upper)) then
        upper_ = upper
      else
        upper_ = 1.0_dp
      end if

      if (present(gap)) then
        gap_ = gap
      else
        gap_ = (1.0_dp,0.0_dp)
      end if
    end subroutine
  end subroutine

  subroutine dos_bilayer()
    continue
  end subroutine

  subroutine dos_trilayer()
    continue
  end subroutine

  subroutine dos_multilayer()
    continue
  end subroutine


  subroutine critical_bilayer(s, m, iterations, stabilization)
    class(superconductor), intent(inout) :: s
    class(conductor),      intent(inout) :: m
    integer,               intent(in   ) :: iterations
    integer,               intent(in   ) :: stabilization
  
    real(dp)                             :: gap   = 0.001_dp
    real(dp)                             :: lower = 0.000_dp
    real(dp)                             :: upper = 1.000_dp
    integer                              :: n
  
    do n=1,iterations
      ! Reset both materials to a weakly superconducting state
      call s%initialize( gap = cmplx(gap,0.000_dp,kind=dp) )
      call m%initialize( gap = cmplx(gap,0.000_dp,kind=dp) )
  
      ! Set the temperature to the midpoint of the current search space
      s%temperature = (lower + upper)/2.0_dp
  
      ! Update the state of the hybrid system
      !do m=1,stabilization
        call m%update
        call s%update
      !  call update(s)
      !end do
  
      ! Check whether the gap has increased, and update the search space
      if (abs(s%get_gap_mean()) >= gap) then
        lower = s%temperature
      else
        upper = s%temperature
      end if
    end do
  end subroutine

end module
