!> Author:   Jabir Ali Ouassou
!> Category: Foundation
!>
!> This module contains constants and functions for doing low-level numerics
!> in Fortran. This includes standard real and complex precisions (sp, dp, qp),
!> a "working precision" (wp) that is used as the default precision in GENEUS;
!> some common numerical constants (like pi and the machine epsilon); as well
!> as basic functions for e.g. constructing and deconstructing complex numbers.

module basic_m
    use :: iso_fortran_env

    ! Declare floating-point precisions
    integer, parameter :: sp = real32  !! Single precision
    integer, parameter :: dp = real64  !! Double precision
    integer, parameter :: qp = real128 !! Quadruple precision
    integer, parameter :: wp = dp      !! Working precision

    ! Define common mathematical constants
    real(wp), parameter :: inf = huge(1.0_wp)        !! Infinity
    real(wp), parameter :: eps = epsilon(1.0_wp)     !! Infinitesimal
    real(wp), parameter :: pi  = atan(1.0_wp)*4.0_wp !! Circle constant
contains
    elemental function re(z) result(x)
    !!  Returns the real part of a complex number z=x+iy.
    !!
    !!  @NOTE:
    !!    This should be replaced by z%re when compilers support it.
        complex(wp), intent(in) :: z !! Complex number
        real(wp)                :: x !! Real part

        x = real(z, kind=wp)
    end function

    elemental function im(z) result(y)
    !!  Returns the imaginary part of a complex number z=x+iy.
    !!
    !!  @NOTE:
    !!    This should be replaced by z%im when compilers support it.
        complex(wp), intent(in) :: z !! Complex number
        real(wp)                :: y !! Imaginary part

        y = aimag(z)
    end function

    elemental function cx(x, y) result(z)
    !!  Returns the complex number z=x+iy.
    !!
    !!  @NOTE:
    !!    This should be rewritten via z%re, z%im when compilers support it.
        real(wp), intent(in)           :: x !! Real part
        real(wp), intent(in), optional :: y !! Imaginary part
        complex(wp)                    :: z !! Complex number

        if (present(y)) then
            z = cmplx(x, y, kind=wp)
        else
            z = cmplx(x, kind=wp)
        end if
    end function

    elemental function arg(z) result(t)
    !!  Returns the complex argument θ of a complex number z=r·exp(iθ).
    !!
    !!  @NOTE:
    !!    This should be rewritten via z%re, z%im when compilers support it.
        complex(wp), intent(in) :: z !! Complex number
        real(wp)                :: t !! Complex argument

        t = atan2(im(z), re(z))
    end function

    pure function unitvector(v) result(r)
    !!  If the argument has a finite norm, then it will be rescaled to a unit
    !!  vector. If that norm is zero, then a zero vector is returned instead.
        real(wp), dimension(3), intent(in) :: v !! Input vector
        real(wp), dimension(3)             :: r !! Unit vector

        r = v/(norm2(v) + eps)
    end function

    pure function nonzero(v) result(r)
    !!  Checks whether or not the argument has a finite norm.
        real(wp), dimension(:), intent(in) :: v !! Input vector
        logical                            :: r !! Conclusion

        r = norm2(v) > eps
    end function
end module
