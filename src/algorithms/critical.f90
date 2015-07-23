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

    ! Perform the binary search
    do n = 1,bisections_
      ! Initialize a weakly superconducting state
      call initialize_all(material, gap_)
      call material%set_gap(gap_)

      ! Update the state of the multilayer structure
      do m = 1,iterations_
        call update_all(material)
        print *,n,m
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
  end subroutine

  !subroutine critical_bilayer(s, m, iterations, stabilization)
  !  class(superconductor), intent(inout) :: s
  !  class(conductor),      intent(inout) :: m
  !  integer,               intent(in   ) :: iterations
  !  integer,               intent(in   ) :: stabilization
  !
  !  real(dp)                             :: gap   = 0.001_dp
  !  real(dp)                             :: lower = 0.000_dp
  !  real(dp)                             :: upper = 1.000_dp
  !  integer                              :: n
  !
  !  do n=1,iterations
  !    ! Reset both materials to a weakly superconducting state
  !    call s%initialize( gap = cmplx(gap,0.000_dp,kind=dp) )
  !    call m%initialize( gap = cmplx(gap,0.000_dp,kind=dp) )
  !
  !    ! Set the temperature to the midpoint of the current search space
  !    s%temperature = (lower + upper)/2.0_dp
  !
  !    ! Update the state of the hybrid system
  !    !do m=1,stabilization
  !      call m%update
  !      call s%update
  !    !  call update(s)
  !    !end do
  !
  !    ! Check whether the gap has increased, and update the search space
  !    if (abs(s%get_gap_mean()) >= gap) then
  !      lower = s%temperature
  !    else
  !      upper = s%temperature
  !    end if
  !  end do
  !end subroutine

end module
