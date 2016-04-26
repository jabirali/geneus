!> Author:   Jabir Ali Ouassou
!> Date:     2015-09-23
!> Category: Foundation
!>
!> This file defines a module containing the machine size of single-precision, double-precision, and quadruple-precision
!> floating point numbers; to declare the floating point precision of a variable, use real(sp), real(dp), or real(qp) as
!> the type of the variable. It also defines the working-precision, which will be the default kind for module procedures.
!> As for module procedures, this library defines some common utility functions for working with e.g. complex numbers.

module math_m
  use :: iso_fortran_env
  use :: stdio_m
  private

  ! Declare which routines to export
  public :: unitvector, re, im, cx

  ! Declare floating-point precisions
  integer,  parameter, public :: sp  = real32              !! Single precision
  integer,  parameter, public :: dp  = real64              !! Double precision
  integer,  parameter, public :: qp  = real128             !! Quadruple precision
  integer,  parameter, public :: wp  = dp                  !! Working precision

  ! Define common mathematical constants
  real(wp), parameter, public :: inf = huge(1.0_wp)        !! Numerical infinity
  real(wp), parameter, public :: eps = epsilon(1.0_wp)     !! Numerical infinitesimal
  real(wp), parameter, public :: pi  = atan(1.0_wp)*4.0_wp !! Circle constant
contains
  pure elemental function re(z) result(x)
    !! Returns the real part of a complex number z=x+iy.
    !! @NOTE This function should be replaced by the structure notation z%re when common compilers support it.
    complex(wp), intent(in) :: z   !! Complex number
    real(wp)                :: x   !! Real part

    x = real(z,kind=wp)
  end function

  pure elemental function im(z) result(y)
    !! Returns the imaginary part of a complex number z=x+iy.
    !! @NOTE This function should be replaced by the structure notation z%im when common compilers support it.
    complex(wp), intent(in) :: z   !! Complex number
    real(wp)                :: y   !! Imaginary part

    y = aimag(z)
  end function

  pure elemental function cx(x,y) result(z)
    !! Returns the complex number z=x+iy.
    real(wp), intent(in)           :: x   !! Real part
    real(wp), intent(in), optional :: y   !! Imaginary part
    complex(wp)                    :: z   !! Complex number

    if (present(y)) then
      z = cmplx(x,y,kind=wp)
    else
      z = cmplx(x,0,kind=wp)
    end if
  end function

  pure function unitvector(v) result(r)
    !! If the argument has a finite norm, then it will be rescaled to a unit
    !! vector. If that norm is zero, then a zero vector is returned instead.
    real(wp), intent(in) :: v(3)   !! Vector
    real(wp)             :: r(3)   !! Unit vector

    r = v/(norm2(v)+eps)
  end function
end module
