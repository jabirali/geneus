! This module defines assert functions, which are used to test the data structures.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-13
! Updated: 2015-07-13


module module_assert
  use module_precision
  use module_spin
  implicit none

  interface assert
    module procedure assert_real, assert_spin
  end interface

  contains
    subroutine assert_print_success(msg)
      character(*),  intent(in) :: msg
      character(50)             :: str

      str = msg
      print *,str,achar(27),'[32;1mSUCCESS',achar(27),'[0m'
    end subroutine

    subroutine assert_print_failure(msg)
      character(*),  intent(in) :: msg
      character(50)             :: str

      str = msg
      print *,str,achar(27),'[31;1mFAILURE',achar(27),'[0m'
    end subroutine

    subroutine assert_real(expr, msg)
      real(dp),      intent(in) :: expr
      character(*),  intent(in) :: msg

      if (expr < 1e-8) then
        call assert_print_failure(msg)
      else
        call assert_print_success(msg)
      end if
    end subroutine

    subroutine assert_spin(a,b,msg)
      type(spin),    intent(in) :: a, b
      character(*),  intent(in) :: msg

      call assert(maxval(abs(a%matrix-b%matrix)), msg)
    end subroutine
end module
