!> Author:   Jabir Ali Ouassou
!> Date:     2015-04-26
!> Category: Foundation
!>
!> This file defines functions that perform some common calculus operations.

module calculus_m
  use :: math_m
  private

  ! Declare which routines to export
  public :: mean, differentiate, integrate, interpolate, linspace

  interface linspace
    !! Public interface for routines that initialize arrays that span certain number ranges.
    module procedure linspace_array_re
  end interface

  interface mean
    !! Public interface for routines that calculate the mean value.
    module procedure mean_array_re, mean_array_cx 
  end interface

  interface differentiate
    !! Public interface for various differentiation routines.
    module procedure differentiate_array_re, differentiate_array_cx
  end interface

  ! Declare public interfaces
  interface integrate
    !! Public interface for various integration routines.
    module procedure integrate_array_re, integrate_array_cx, &
                     integrate_range_re, integrate_range_cx
  end interface

  interface interpolate
    !! Public interface for various interpolation routines.
    module procedure interpolate_array_re, interpolate_array_cx, &
                     interpolate_point_re, interpolate_point_cx
  end interface 
contains

  !--------------------------------------------------------------------------------
  ! Specific implementations of the `linspace` interface
  !--------------------------------------------------------------------------------

  pure subroutine linspace_array_re(array, first, last)
    !! Populates an existing array with elements from `first` to `last`, inclusive.
    real(wp), dimension(:), intent(out) :: array   !! Output array to populate
    real(wp),               intent(in)  :: first   !! Value of first element
    real(wp),               intent(in)  :: last    !! Value of last  element
    integer                             :: n

    do n=1,size(array)
      array(n) = first + ((last-first)*(n-1))/(size(array)-1)
    end do
  end subroutine

  !--------------------------------------------------------------------------------
  ! Specific implementations of the `mean` interface
  !--------------------------------------------------------------------------------

  pure function mean_array_re(x) result(r)
    !! Calculates the mean value of a real-valued array.
    real(wp), dimension(:), intent(in) :: x   !! Real-valued array
    real(wp)                           :: r   !! Mean value <x>

    r = sum(x)/max(1,size(x))
  end function

  pure function mean_array_cx(x) result(r)
    !! Calculates the mean value of a complex-valued array.
    complex(wp), dimension(:), intent(in) :: x  !! Complex-valued array
    complex(wp)                           :: r  !! Mean value <x>

    r = sum(x)/max(1,size(x))
  end function

  !--------------------------------------------------------------------------------
  ! Specific implementations of the `differentiate` interface
  !--------------------------------------------------------------------------------

  pure function differentiate_array_re(x, y) result(r)
    !! This function calculates the numerical derivative of an array y with respect to x, using a central difference approximation
    !! at the interior points and forward/backward difference approximations at the exterior points. Note that since all the three
    !! approaches yield two-point approximations of the derivative, the mesh spacing of x does not necessarily have to be uniform.
    real(wp), dimension(:),       intent(in)  :: x   !! Variable x
    real(wp), dimension(size(x)), intent(in)  :: y   !! Function y(x)
    real(wp), dimension(size(x))              :: r   !! Derivative dy/dx

    ! Differentiate using finite differences
    associate(n => size(x))
      r(   1   ) = (y( 1+1 ) - y(  1  ))/(x( 1+1 ) - x(  1  ))
      r(1+1:n-1) = (y(1+2:n) - y(1:n-2))/(x(1+2:n) - x(1:n-2))
      r(   n   ) = (y(  n  ) - y( n-1 ))/(x(  n  ) - x( n-1 ))
    end associate
  end function

  pure function differentiate_array_cx(x, y) result(r)
    !! Complex version of differentiate_array_re.
    real(wp),    dimension(:),       intent(in)  :: x   !! Variable x
    complex(wp), dimension(size(x)), intent(in)  :: y   !! Function y(x)
    complex(wp), dimension(size(x))              :: r   !! Derivative dy/dx

    ! Differentiate using finite differences
    associate(n => size(x))
      r(   1   ) = (y( 1+1 ) - y(  1  ))/(x( 1+1 ) - x(  1  ))
      r(1+1:n-1) = (y(1+2:n) - y(1:n-2))/(x(1+2:n) - x(1:n-2))
      r(   n   ) = (y(  n  ) - y( n-1 ))/(x(  n  ) - x( n-1 ))
    end associate
  end function

  !--------------------------------------------------------------------------------
  ! Specific implementations of the `integrate` interface
  !--------------------------------------------------------------------------------

  pure function integrate_array_re(x, y) result(r)
    !! This function calculates the integral of an array y with respect to x using a trapezoid
    !! approximation. Note that the mesh spacing of x does not necessarily have to be uniform.
    real(wp), dimension(:),       intent(in)  :: x   !! Variable x
    real(wp), dimension(size(x)), intent(in)  :: y   !! Function y(x)
    real(wp)                                  :: r   !! Integral ∫y(x)·dx

    ! Integrate using the trapezoidal rule
    associate(n => size(x))
      r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
    end associate
  end function

  pure function integrate_array_cx(x, y) result(r)
    !! Complex version of integrate_array_re.
    real(wp),    dimension(:),       intent(in) :: x   !! Variable x
    complex(wp), dimension(size(x)), intent(in) :: y   !! Function y(x)
    complex(wp)                                 :: r   !! Integral ∫y(x)·dx

    ! Integrate using the trapezoidal rule
    associate(n => size(x))
      r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
    end associate
  end function

  function integrate_range_re(x, y, a, b) result(r)
    !! This function constructs a piecewise hermitian cubic interpolation of an array y(x) based on
    !! discrete numerical data, and subsequently evaluates the integral of the interpolation in the
    !! range (a,b). Note that the mesh spacing of x does not necessarily have to be uniform.
    real(wp), dimension(:),       intent(in)  :: x  !! Variable x
    real(wp), dimension(size(x)), intent(in)  :: y  !! Function y(x)
    real(wp),                     intent(in)  :: a  !! Left endpoint
    real(wp),                     intent(in)  :: b  !! Right endpoint
    real(wp)                                  :: r  !! Integral ∫y(x)·dx

    external                                  :: dpchez
    real(wp), external                        :: dpchqa
    real(wp), dimension(size(x))              :: d
    integer                                   :: err

    ! Create a PCHIP interpolation of the input data
    call dpchez(size(x), x, y, d, .false., 0, 0, err)

    ! Integrate the interpolation in the provided range
    r = dpchqa(size(x), x, y, d, a, b, err)
  end function

  function integrate_range_cx(x, y, a, b) result(r)
    !! Complex version of integrate_range_re.
    real(wp),    dimension(:),       intent(in)  :: x   !! Variable x
    complex(wp), dimension(size(x)), intent(in)  :: y   !! Function y(x)
    real(wp),                        intent(in)  :: a   !! Left endpoint
    real(wp),                        intent(in)  :: b   !! Right endpoint
    complex(wp)                                  :: r   !! Integral ∫y(x)·dx

    ! Integrate the real and imaginary parts separately
    r = cx( integrate_range_re(x, re(y), a, b),&
            integrate_range_re(x, im(y), a, b) )
  end function

  !--------------------------------------------------------------------------------
  ! Specific implementations of the `interpolate` interface
  !--------------------------------------------------------------------------------

  function interpolate_array_re(x, y, p) result(r)
    !! This function constructs a piecewise hermitian cubic interpolation of an array y(x) based on discrete numerical data,
    !! and evaluates the interpolation at points p. Note that the mesh spacing of x does not necessarily have to be uniform.
    real(wp), dimension(:),       intent(in)  :: x   !! Variable x
    real(wp), dimension(size(x)), intent(in)  :: y   !! Function y(x)
    real(wp), dimension(:),       intent(in)  :: p   !! Interpolation domain p
    real(wp), dimension(size(p))              :: r   !! Interpolation result y(p)

    external                                  :: dpchez
    external                                  :: dpchfe
    real(wp), dimension(size(x))              :: d
    integer                                   :: err

    ! Create a PCHIP interpolation of the input data
    call dpchez(size(x), x, y, d, .false., 0, 0, err)

    ! Extract the interpolated data at provided points
    call dpchfe(size(x), x, y, d, 1, .false., size(p), p, r, err)
  end function

  function interpolate_array_cx(x, y, p) result(r)
    !! Complex version of interpolate_array_re.
    real(wp),    dimension(:),       intent(in)  :: x   !! Variable x
    complex(wp), dimension(size(x)), intent(in)  :: y   !! Function y(x)
    real(wp),    dimension(:),       intent(in)  :: p   !! Interpolation domain p
    complex(wp), dimension(size(p))              :: r   !! Interpolation result y(p)

    ! Interpolate the real and imaginary parts separately
    r = cx( interpolate_array_re(x, re(y), p),&
            interpolate_array_re(x, im(y), p) )
  end function

  function interpolate_point_re(x, y, p) result(r)
    !! Wrapper for interpolate_pchip that accepts scalar arguments.
    real(wp), dimension(:),       intent(in)  :: x   !! Variable x
    real(wp), dimension(size(x)), intent(in)  :: y   !! Function y(x)
    real(wp),                     intent(in)  :: p   !! Interpolation point p
    real(wp)                                  :: r   !! Interpolation result y(p)
    real(wp)                                  :: rs(1) 

    ! Perform the interpolation
    rs = interpolate_array_re(x, y, [p])

    ! Extract the scalar result
    r = rs(1)
  end function

  function interpolate_point_cx(x, y, p) result(r)
    !! Complex version of interpolate_point_re.
    real(wp),    dimension(:),       intent(in)  :: x   !! Variable x
    complex(wp), dimension(size(x)), intent(in)  :: y   !! Function y(x)
    real(wp),                        intent(in)  :: p   !! Interpolation point p
    complex(wp)                                  :: r   !! Interpolation result y(p)
    complex(wp)                                  :: rs(1) 

    ! Perform the interpolation
    rs = cx( interpolate_array_re(x, re(y), [p]),&
             interpolate_array_re(x, im(y), [p]) )

    ! Extract the scalar result
    r = rs(1)
  end function
end module
