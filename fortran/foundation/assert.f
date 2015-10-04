! This module defines the interfaces 'check' and 'assert', which are used for unit testing of other modules.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-13
! Updated: 2015-10-04


module mod_assert
  use mod_system
  use mod_math, only: eps, wp
  use mod_spin, only: spin
  implicit none

  interface check
    module procedure check_scalar, check_array
  end interface

  interface assert
    module procedure assert_scalar, assert_array
  end interface
contains
  impure subroutine assert_print(valid, msg)
    ! Prints to standard out whether an assertion succeeded or failed.
    logical,      intent(in) :: valid
    character(*), intent(in) :: msg
    character(60)            :: str

    str = msg
    if (valid) then
      print *,' :: ', str, color_green, 'SUCCESS', color_none
    else
      print *,' :: ', str, color_red,   'FAILURE', color_none
    end if
  end subroutine

  impure subroutine assert_scalar(expr, msg)
    ! Implements an assert function for polymorphic scalar arguments.
    class(*),     intent(in) :: expr
    character(*), intent(in) :: msg

    call assert_print(check(expr), msg)
  end subroutine

  impure subroutine assert_array(expr, msg)
    ! Implements an assert function for polymorphic array arguments.
    class(*),     intent(in) :: expr(:)
    character(*), intent(in) :: msg

    call assert_print(check(expr), msg)
  end subroutine

  pure function check_array(expr) result(r)
    ! Implements a polymorphic comparison function for array arguments, which returns .true. when all expressions are numerically zero.
    class(*), intent(in) :: expr(:)
    integer              :: n
    logical              :: r

    do n=1,size(expr)
      r = check_scalar(expr(n))
      if (.not. r) return
    end do
  end function

  pure recursive function check_scalar(expr) result(r)
    ! Implements a polymorphic comparison function for scalar arguments, which returns .true. when the expression is numerically zero.
    class(*), intent(in) :: expr
    logical              :: r

    select type(expr)
      type is (logical)
        r = expr
      type is (integer)
        r = check(expr == 0)
      type is (real(wp))
        r = check(abs(expr) < sqrt(eps))
      type is (complex(wp))
        r = check(abs(expr) < sqrt(eps))
      type is(spin)
        r = check(reshape(expr%matrix,[4]))
    end select
  end function
end module
