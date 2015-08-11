! This module defines the interface 'assert', which is used for unit testing.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-13
! Updated: 2015-08-09


module mod_assert
  use mod_system
  use mod_spin
  implicit none

  interface assert
    module procedure assert_logical, assert_integer, assert_real, assert_complex, assert_spin, &
                     assert_logical_array, assert_integer_array, assert_real_array, assert_complex_array
  end interface

  contains
    subroutine assert_logical(expr, msg)
      logical,      intent(in) :: expr
      character(*), intent(in) :: msg
      character(60)            :: str

      str = msg
      if (colors) then
        if (expr) then
          print *,' :: ', str, color_green, 'SUCCESS', color_none
        else
          print *,' :: ', str, color_red,   'FAILURE', color_none
        end if
      else
        if (expr) then
          print *,' :: ', str, 'SUCCESS'
        else
          print *,' :: ', str, 'FAILURE'
        end if
      end if
    end subroutine

    subroutine assert_logical_array(expr, msg)
      logical,      intent(in) :: expr(:)
      character(*), intent(in) :: msg

      call assert(all(expr), msg)
    end subroutine

    subroutine assert_integer(expr, msg)
      integer,      intent(in) :: expr
      character(*), intent(in) :: msg

      call assert(expr == 0, msg)
    end subroutine

    subroutine assert_integer_array(expr, msg)
      integer,      intent(in) :: expr(:)
      character(*), intent(in) :: msg

      call assert(all(expr == 0), msg)
    end subroutine

    subroutine assert_real(expr, msg)
      real(dp),     intent(in) :: expr
      character(*), intent(in) :: msg

      call assert(abs(expr) < 1e-8, msg)
    end subroutine

    subroutine assert_real_array(expr, msg)
      real(dp),     intent(in) :: expr(:)
      character(*), intent(in) :: msg

      call assert(all(abs(expr) < 1e-8), msg)
    end subroutine

    subroutine assert_complex(expr, msg)
      complex(dp),  intent(in) :: expr
      character(*), intent(in) :: msg

      call assert(abs(expr) < 1e-8, msg)
    end subroutine

    subroutine assert_complex_array(expr, msg)
      complex(dp),  intent(in) :: expr(:)
      character(*), intent(in) :: msg

      call assert(all(abs(expr) < 1e-8), msg)
    end subroutine

    subroutine assert_spin(expr,msg)
      type(spin),    intent(in) :: expr
      character(*),  intent(in) :: msg

      call assert(maxval(abs(expr%matrix)), msg)
    end subroutine
end module
