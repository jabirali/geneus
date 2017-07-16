!> Author:   Jabir Ali Ouassou
!> Category: Foundation
!>
!> This file defines a module containing the machine size of single-precision, double-precision, and quadruple-precision
!> floating point numbers; to declare the floating point precision of a variable, use real(sp), real(dp), or real(qp) as
!> the type of the variable. It also defines the working-precision, which will be the default kind for module procedures.
!> As for module procedures, this library defines some common utility functions for working with e.g. complex numbers.
!>
!> @TODO:
!>   To make future reuse of the math routines in other projects easier, and also clean up the s4tran code 
!>   internally, the math routines should be refactored as follows:
!>    (1) This module should be renamed to e.g. flop_m, because it has to do with basic floating-point operations;
!>    (2) All the other math modules matrix_m, calculus_m, etc. should then import flop_m and do their thing;
!>    (3) A new module with the name math_m imports all these other modules, and reexports their function interfaces;
!>    (4) If and when parametrized derived types get implemented in GFortran too, spin_m and nambu_m can be merged into
!>        a general matrix library for N×N stuff. Until then, they stay out of the general math library I'm making.
!>   This way, other libraries and programs can simply "use :: math_m", and automatically get access to everything,
!>   without having to cherry-pick the parts of the library that are relevant to import. (All of math_m is linked
!>   together into one library libs4tran.a anyway, so this shouldn't affect compilation time. The only thing that
!>   would speed up compilation time would be to adopt submodules, but that implies dropping gfortran 5.x support,
!>   and at the same time adopting the ugly C habit of defining the interface in the .h file away from the function.)

module math_m
  use :: iso_fortran_env
  use :: stdio_m
  private

  ! Declare which routines to export
  public :: unitvector, re, im, cx, arg

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

  pure elemental function arg(z) result(t)
    !! Returns the complex argument θ of a complex number z=r·exp(iθ). 
    !! @NOTE: this function returns the normalized angle θ/π.
    complex(wp), intent(in) :: z   !! Complex number
    real(wp)                :: t   !! Complex argument

    t = atan2(im(z), re(z))/pi
  end function

  pure function unitvector(v) result(r)
    !! If the argument has a finite norm, then it will be rescaled to a unit
    !! vector. If that norm is zero, then a zero vector is returned instead.
    real(wp), intent(in) :: v(3)   !! Vector
    real(wp)             :: r(3)   !! Unit vector

    r = v/(norm2(v)+eps)
  end function
end module
