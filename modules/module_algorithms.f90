!

module module_algorithms
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
      do m=1,stabilization
        call m%update
        call s%update
      end do
  
      ! Check whether the gap has increased, and update the search space
      if (abs(s%get_gap_mean()) >= gap) then
        lower = s%temperature
      else
        upper = s%temperature
      end if
    end do
  end subroutine
end module
