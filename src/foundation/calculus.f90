!> Author:   Jabir Ali Ouassou
!> Category: Foundation
!>
!> This module defines some functions that perform common calculus operations,
!> including differentiating, integrating, and interpolating discretized data.

module calculus_m
    use :: basic_m
contains
    !---------------------------------------------------------------------------
    ! Specific implementations of the `mean` interface
    !---------------------------------------------------------------------------

    pure function mean_array_re(x) result(r)
    !!  Calculates the mean value of a real-valued array.
        real(wp), dimension(:), intent(in) :: x !! Real-valued array
        real(wp)                           :: r !! Mean value <x>

        r = sum(x)/max(1, size(x))
    end function

    pure function mean_array_cx(x) result(r)
    !!  Calculates the mean value of a complex-valued array.
        complex(wp), dimension(:), intent(in) :: x !! Complex-valued array
        complex(wp)                           :: r !! Mean value <x>

        r = sum(x)/max(1, size(x))
    end function

    !---------------------------------------------------------------------------
    ! Specific implementations of the `differentiate` interface
    !---------------------------------------------------------------------------

    pure function differentiate_array_re(x, y) result(r)
    !!  Calculates the numerical derivative of an array y wrt. x using central
    !!  differences at the interior points and forward/backward differences at
    !!  the exterior points. All three approaches yield two-point approximations
    !!  of the derivatives, thus the mesh spacing does not have to be uniform.
        real(wp), dimension(:),       intent(in) :: x !! Variable x
        real(wp), dimension(size(x)), intent(in) :: y !! Function y(x)
        real(wp), dimension(size(x))             :: r !! Derivative dy/dx

        ! Differentiate using finite differences
        associate (n => size(x))
            r(1)     = (y(2) - y(1)) / (x(2) - x(1))
            r(2:n-1) = (y(3:n) - y(1:n-2)) / (x(3:n) - x(1:n-2))
            r(n)     = (y(n) - y(n-1)) / (x(n) - x(n-1))
        end associate
    end function

    pure function differentiate_array_cx(x, y) result(r)
    !!  Complex version of differentiate_array_re.
        real(wp),    dimension(:),       intent(in) :: x !! Variable x
        complex(wp), dimension(size(x)), intent(in) :: y !! Function y(x)
        complex(wp), dimension(size(x))             :: r !! Derivative dy/dx

        ! Differentiate using finite differences
        associate (n => size(x))
            r(1)     = (y(2) - y(1)) / (x(2) - x(1))
            r(2:n-1) = (y(3:n) - y(1:n-2)) / (x(3:n) - x(1:n-2))
            r(n)     = (y(n) - y(n-1)) / (x(n) - x(n-1))
        end associate
    end function

    !---------------------------------------------------------------------------
    ! Specific implementations of the `integrate` interface
    !---------------------------------------------------------------------------

    pure function integrate_array_re(x, y) result(r)
    !!  Calculates the integral of an array y wrt. x using the trapezoid
    !!  method. The mesh spacing does not necessarily have to be uniform.
        real(wp), dimension(:),       intent(in) :: x !! Variable x
        real(wp), dimension(size(x)), intent(in) :: y !! Function y(x)
        real(wp)                                 :: r !! Integral ∫y(x)·dx

        ! Integrate using the trapezoidal rule
        associate (n => size(x))
            r = sum((y(2:n) + y(1:n-1)) * (x(2:n) - x(1:n-1)))/2
        end associate
    end function

    pure function integrate_array_cx(x, y) result(r)
    !!  Complex version of integrate_array_re.
        real(wp),    dimension(:),       intent(in) :: x !! Variable x
        complex(wp), dimension(size(x)), intent(in) :: y !! Function y(x)
        complex(wp)                                 :: r !! Integral ∫y(x)·dx

        ! Integrate using the trapezoidal rule
        associate (n => size(x))
            r = sum((y(2:n) + y(1:n-1)) * (x(2:n) - x(1:n-1)))/2
        end associate
    end function

    function integrate_range_re(x, y, a, b) result(r)
    !!  Constructs a piecewise hermitian cubic interpolation of an array y(x)
    !!  from discrete numerical data, and then integrates the interpolation in
    !!  the range (a, b). The mesh spacing does not have to be uniform.
        real(wp), dimension(:),       intent(in) :: x !! Variable x
        real(wp), dimension(size(x)), intent(in) :: y !! Function y(x)
        real(wp), intent(in)                     :: a !! Left endpoint
        real(wp), intent(in)                     :: b !! Right endpoint
        real(wp)                                 :: r !! Integral ∫y(x)·dx

        external                     :: dpchez
        real(wp), external           :: dpchqa
        real(wp), dimension(size(x)) :: d
        integer                      :: err

        ! Create a PCHIP interpolation of the input data
        call dpchez(size(x), x, y, d, .false., 0, 0, err)

        ! Integrate the interpolation in the provided range
        r = dpchqa(size(x), x, y, d, a, b, err)
    end function

    function integrate_range_cx(x, y, a, b) result(r)
    !!  Complex version of integrate_range_re.
        real(wp),    dimension(:),       intent(in) :: x !! Variable x
        complex(wp), dimension(size(x)), intent(in) :: y !! Function y(x)
        real(wp),    intent(in)                     :: a !! Left endpoint
        real(wp),    intent(in)                     :: b !! Right endpoint
        complex(wp)                                 :: r !! Integral ∫y(x)·dx

        ! Integrate the real and imaginary parts separately
        ! TODO: Replace this with two sequential integrations.
        r = cx(integrate_range_re(x, re(y), a, b), &
               integrate_range_re(x, im(y), a, b))
    end function

    !---------------------------------------------------------------------------
    ! Specific implementations of the `interpolate` interface
    !---------------------------------------------------------------------------

    function interpolate_array_re(x, y, p) result(r)
    !!  Constructs a piecewise hermitian cubic interpolation of an array y(x)
    !!  based on discrete numerical data and evaluates the interpolation at p.
    !!  Note that the mesh spacing does not necessarily have to be uniform.
        real(wp), dimension(:),       intent(in) :: x !! Variable x
        real(wp), dimension(size(x)), intent(in) :: y !! Function y(x)
        real(wp), dimension(:),       intent(in) :: p !! Point array p
        real(wp), dimension(size(p))             :: r !! Interpolation y(p)

        external                     :: dpchez
        external                     :: dpchfe
        real(wp), dimension(size(x)) :: d
        integer                      :: err

        ! Create a PCHIP interpolation of the input data
        call dpchez(size(x), x, y, d, .false., 0, 0, err)

        ! Extract the interpolated data at provided points
        call dpchfe(size(x), x, y, d, 1, .false., size(p), p, r, err)
    end function

    function interpolate_array_cx(x, y, p) result(r)
    !!  Complex version of interpolate_array_re.
        real(wp),    dimension(:),       intent(in) :: x !! Variable x
        complex(wp), dimension(size(x)), intent(in) :: y !! Function y(x)
        real(wp),    dimension(:),       intent(in) :: p !! Point array p
        complex(wp), dimension(size(p))             :: r !! Interpolation y(p)

        ! Interpolate the real and imaginary parts separately
        ! TODO: Replace this with two sequential interpolations.
        r = cx(interpolate_array_re(x, re(y), p), &
               interpolate_array_re(x, im(y), p))
    end function

    function interpolate_point_re(x, y, p) result(r)
    !!  Wrapper for interpolate_array_re that accepts scalar arguments.
        real(wp), dimension(:),       intent(in) :: x !! Variable x
        real(wp), dimension(size(x)), intent(in) :: y !! Function y(x)
        real(wp), intent(in)                     :: p !! Single point p
        real(wp)                                 :: r !! Interpolation y(p)

        real(wp), dimension(1) :: rs

        ! Perform the interpolation
        rs = interpolate_array_re(x, y, [p])

        ! Extract the scalar result
        r = rs(1)
    end function

    function interpolate_point_cx(x, y, p) result(r)
    !!  Complex version of interpolate_point_re.
        real(wp),    dimension(:),       intent(in) :: x !! Variable x
        complex(wp), dimension(size(x)), intent(in) :: y !! Function y(x)
        real(wp),                        intent(in) :: p !! Single point p
        complex(wp)                                 :: r !! Interpolation y(p)

        complex(wp), dimension(1) :: rs

        ! Perform the interpolation
        ! TODO: Replace this with two sequential interpolations.
        rs = cx(interpolate_array_re(x, re(y), [p]), &
                interpolate_array_re(x, im(y), [p]))

        ! Extract the scalar result
        r = rs(1)
    end function

    pure function interpolate_point_matrix_re(x, y, p) result(r)
    !!  Interpolates a matrix function using Catmull-Rom splines.
        real(wp), dimension(:),          intent(in) :: x !! Variable x
        real(wp), dimension(:, :, :),    intent(in) :: y !! Function y(x)
        real(wp),                        intent(in) :: p !! Single point p
        real(wp), dimension(size(y, 1), size(y, 2)) :: r !! Interpolation y(p)

        integer  :: n, m
        real(wp) :: t

        ! Find the nearest known point y(x)
        m = size(x)
        n = floor(p*(m - 1) + 1); 

        ! Perform the interpolation
        if (n <= 0) then
            ! Exterior: nearest extrapolation
            r = y(:, :, 1)
        else if (n >= m) then
            ! Exterior: nearest extrapolation
            r = y(:, :, m)
        else
            ! Interior: spline interpolation
            t = (p - x(max(n, 1)))/(x(min(n + 1, m)) - x(max(n, 1)))
            r = y(:, :, max(n - 1, 1)) * (     -0.5*t +1.0*t**2 -0.5*t**3) &
              + y(:, :, max(n - 0, 1)) * (+1.0        -2.5*t**2 +1.5*t**3) &
              + y(:, :, min(n + 1, m)) * (     +0.5*t +2.0*t**2 -1.5*t**3) &
              + y(:, :, min(n + 2, m)) * (            -0.5*t**2 +0.5*t**3)
        end if
    end function

    !---------------------------------------------------------------------------
    ! Specific implementations of the `linspace` interface
    !---------------------------------------------------------------------------

    pure subroutine linspace_array_re(array, first, last)
    !!  Populates an array with elements from `first` to `last`, inclusive.
        real(wp), dimension(:), intent(out) :: array !! Output array to populate
        real(wp),               intent(in)  :: first !! Value of first element
        real(wp),               intent(in)  :: last  !! Value of last  element

        integer :: n

        do n = 1, size(array)
            array(n) = first + ((last - first)*(n - 1))/(size(array) - 1)
        end do
    end subroutine
end module
