! This file defines a module containing the machine size of single-precision, double-precision, and quadruple-precision
! floating point numbers; to declare the floating point precision of a variable, use real(sp), real(dp), or real(qp) as
! the type of the variable. It also defines the working-precision, which will be the default kind for module procedures.
! As for module procedures, this library defines common utility functions for working with complex numbers and matrices.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-09-23
! Updated: 2016-03-04

module mod_math
  use :: iso_fortran_env
  use :: mod_stdio
  implicit none
  public

  ! Declare floating-point precisions
  integer,     parameter :: sp  = REAL32              ! Single precision
  integer,     parameter :: dp  = REAL64              ! Double precision
  integer,     parameter :: qp  = REAL128             ! Quadruple precision
  integer,     parameter :: wp  = dp                  ! Working precision

  ! Define common mathematical constants
  real(wp),    parameter :: inf = huge(1.0_wp)        ! Numerical infinity
  real(wp),    parameter :: eps = epsilon(1.0_wp)     ! Numerical infinitesimal
  real(wp),    parameter :: pi  = atan(1.0_wp)*4.0_wp ! Circle constant

  ! Define common identity matrices
  real(wp),    parameter :: mateye2(2,2) = reshape([1,0,0,1],[2,2])
  real(wp),    parameter :: mateye3(3,3) = reshape([1,0,0,0,1,0,0,0,1],[3,3])
  real(wp),    parameter :: mateye4(4,4) = reshape([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],[4,4])

  ! Interfaces for math routines
  interface differentiate
    module procedure differentiate_linear, differentiate_linear_cx
  end interface

  interface integrate
    module procedure integrate_linear, integrate_linear_cx, &
                     integrate_pchip,  integrate_pchip_cx
  end interface

  interface interpolate
    module procedure interpolate_pchip, interpolate_pchip_cx
  end interface
contains

  !---------------------------------------------------------------------------!
  !                        COMPLEX NUMBER PROCEDURES                          !
  !---------------------------------------------------------------------------!

  pure elemental function re(z) result(x)
    ! Returns the real part of a complex number z=x+iy.
    complex(wp), intent(in) :: z
    real(wp)                :: x

    x = real(z,kind=wp)
  end function

  pure elemental function im(z) result(y)
    ! Returns the imaginary part of a complex number z=x+iy.
    complex(wp), intent(in) :: z
    real(wp)                :: y

    y = aimag(z)
  end function

  pure elemental function cx(x,y) result(z)
    ! Returns the complex number z=x+iy.
    real(wp), intent(in) :: x
    real(wp), intent(in) :: y
    complex(wp)          :: z

    z = cmplx(x,y,kind=wp)
  end function

  !---------------------------------------------------------------------------!
  !                       ELEMENTARY MATRIX PROCEDURES                        !
  !---------------------------------------------------------------------------!

  pure subroutine linspace(array, first, last)
    ! Populates an array with elements from 'first' to 'last', inclusive.
    real(wp), intent(inout) :: array(:)
    real(wp), intent(in   ) :: first, last
    integer                 :: n

    do n=1,size(array)
      array(n) = first + ((last-first)*(n-1))/(size(array)-1)
    end do
  end subroutine

  pure function mateye(n) result(A)
    ! Constructs an n×n identity matrix.
    integer, intent(in) :: n
    integer             :: i, j
    real(wp)            :: A(n,n)

    ! Initialize by exploiting integer arithmetic to avoid multiple passes
    do i=1,n
      do j=1,n
        A(j,i) = (i/j)*(j/i)
      end do
    end do
  end function

  impure function matrnd(n) result(C)
    ! Constructs an n×n random matrix.
    integer, allocatable :: seed(:)
    integer              :: n, m, u
    complex(wp)          :: C(n,n)
    real(wp)             :: B(n,n)
    real(wp)             :: A(n,n)

    ! Check the size of a random seed
    call random_seed(size = m)

    ! Allocate memory for the random seed
    allocate(seed(m))

    ! Initialize the seed using random numbers from the operating system
    open(newunit=u, file="/dev/urandom", access="stream", form="unformatted", action="read", status="old")
    read(unit=u) seed
    close(unit=u)

    ! Initialize the random number generator
    call random_seed(put = seed)

    ! Construct a random complex matrix
    call random_number(A)
    call random_number(B)
    C = cx(A,B)

    ! Deallocate dynamic memory
    deallocate(seed)
  end function

  pure function matinv2(A) result(B)
    ! Performs a direct calculation of the inverse of a 2×2 matrix.
    complex(wp), intent(in) :: A(2,2)
    complex(wp)             :: B(2,2)
    complex(wp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)
  end function

  pure function matdivl2(A,B) result(C)
    ! Performs a direct calculation of a 2×2 matrix left division.
    complex(wp), intent(in) :: A(2,2)
    complex(wp), intent(in) :: B(2,2)
    complex(wp)             :: C(2,2)
    complex(wp)             :: detinv

    ! Calculate the inverse determinant of the left matrix
    detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the elements of the resulting matrix
    C(1,1) = detinv * (A(2,2)*B(1,1) - A(1,2)*B(2,1))
    C(2,1) = detinv * (A(1,1)*B(2,1) - A(2,1)*B(1,1))
    C(1,2) = detinv * (A(2,2)*B(1,2) - A(1,2)*B(2,2))
    C(2,2) = detinv * (A(1,1)*B(2,2) - A(2,1)*B(1,2))
  end function

  pure function matdivr2(A,B) result(C)
    ! Performs a direct calculation of a 2×2 matrix right division.
    complex(wp), intent(in) :: A(2,2)
    complex(wp), intent(in) :: B(2,2)
    complex(wp)             :: C(2,2)
    complex(wp)             :: detinv

    ! Calculate the inverse determinant of the right matrix
    detinv = 1/(B(1,1)*B(2,2) - B(1,2)*B(2,1))

    ! Calculate the elements of the resulting matrix
    C(1,1) = detinv * (A(1,1)*B(2,2) - A(1,2)*B(2,1))
    C(2,1) = detinv * (A(2,1)*B(2,2) - A(2,2)*B(2,1))
    C(1,2) = detinv * (A(1,2)*B(1,1) - A(1,1)*B(1,2))
    C(2,2) = detinv * (A(2,2)*B(1,1) - A(2,1)*B(1,2))
  end function

  pure function matinv3(A) result(B)
    ! Performs a direct calculation of the inverse of a 3×3 matrix.
    ! [Based on the subroutine M33INV by David G. Simpson, NASA.]
    complex(wp), intent(in) :: A(3,3)
    complex(wp)             :: B(3,3)
    complex(wp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
  end function

  pure function matinv4(A) result(B)
    ! Performs a direct calculation of the inverse of a 4×4 matrix.
    ! [Based on the subroutine M44INV by David G. Simpson, NASA.]
    complex(wp), intent(in) :: A(4,4)
    complex(wp)             :: B(4,4)
    complex(wp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = &
      1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
       - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
       + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
       - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    ! Calculate the inverse of the matrix
    B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
  end function

  !---------------------------------------------------------------------------!
  !                       ELEMENTARY CALCULUS PROCEDURES                      !
  !---------------------------------------------------------------------------!

  pure function differentiate_linear(x, y) result(r)
    ! This function calculates the numerical derivative of an array y with respect to x, using a central difference approximation
    ! at the interior points and forward/backward difference approximations at the exterior points. Note that since all the three
    ! approaches yield two-point approximations of the derivative, the mesh spacing of x does not necessarily have to be uniform.
    real(wp), intent(in)  :: x(:)
    real(wp), intent(in)  :: y(:)
    real(wp), allocatable :: r(:)
    integer               :: n

    ! Check the size of input arrays
    n = min(size(x), size(y))

    ! Allocate memory for the results
    allocate(r(n))

    ! Differentiate using finite differences
    r(   1   ) = (y( 1+1 ) - y(  1  ))/(x( 1+1 ) - x(  1  ))
    r(1+1:n-1) = (y(1+2:n) - y(1:n-2))/(x(1+2:n) - x(1:n-2))
    r(   n   ) = (y(  n  ) - y( n-1 ))/(x(  n  ) - x( n-1 ))
  end function

  pure function differentiate_linear_cx(x, y) result(r)
    ! Complex version of differentiate_linear.
    real(wp),    intent(in)  :: x(:)
    complex(wp), intent(in)  :: y(:)
    complex(wp), allocatable :: r(:)
    integer                  :: n

    ! Check the size of input arrays
    n = min(size(x), size(y))

    ! Allocate memory for the results
    allocate(r(n))

    ! Differentiate using finite differences
    r(   1   ) = (y( 1+1 ) - y(  1  ))/(x( 1+1 ) - x(  1  ))
    r(1+1:n-1) = (y(1+2:n) - y(1:n-2))/(x(1+2:n) - x(1:n-2))
    r(   n   ) = (y(  n  ) - y( n-1 ))/(x(  n  ) - x( n-1 ))
  end function

  pure function integrate_linear(x, y) result(r)
    ! This function calculates the integral of an array y with respect to x using a trapezoid
    ! approximation. Note that the mesh spacing of x does not necessarily have to be uniform.
    real(wp), intent(in)  :: x(:)
    real(wp), intent(in)  :: y(:)
    real(wp)              :: r
    integer               :: n

    ! Check the size of input arrays
    n = min(size(x), size(y))

    ! Integrate using the trapezoidal rule
    r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
  end function

  pure function integrate_linear_cx(x, y) result(r)
    ! Complex version of integrate_linear.
    real(wp),    intent(in) :: x(:)
    complex(wp), intent(in) :: y(:)
    complex(wp)             :: r
    integer                 :: n

    ! Check the size of input arrays
    n = min(size(x), size(y))

    ! Integrate using the trapezoidal rule
    r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
  end function

  function integrate_pchip(x, y, a, b) result(r)
    ! This function constructs a piecewise hermitian cubic interpolation of an array y(x) based on
    ! discrete numerical data, and subsequently evaluates the integral of the interpolation in the
    ! range (a,b). Note that the mesh spacing of x does not necessarily have to be uniform.
    real(wp), external    :: dpchqa
    real(wp), intent(in)  :: x(:)
    real(wp), intent(in)  :: y(:)
    real(wp), allocatable :: d(:)
    real(wp), intent(in)  :: a
    real(wp), intent(in)  :: b
    real(wp)              :: r
    integer               :: err
    integer               :: n

    ! Check the size of input arrays
    n = min(size(x), size(y))

    ! Allocate memory for the derivatives
    allocate(d(n))

    ! Create a PCHIP interpolation of the input data
    call dpchez(n, x, y, d, .false., 0, 0, err)

    ! Integrate the interpolation in the provided range
    r = dpchqa(n, x, y, d, a, b, err)

    ! Clean up workspace memory
    deallocate(d)
  end function

  function integrate_pchip_cx(x, y, a, b) result(r)
    ! Wrapper for integrate_pchip that accepts complex arguments.
    real(wp),    intent(in)  :: x(:)
    complex(wp), intent(in)  :: y(:)
    real(wp),    intent(in)  :: a
    real(wp),    intent(in)  :: b
    complex(wp)              :: r

    ! Integrate the real and imaginary parts separately
    r = cx( integrate_pchip(x, re(y), a, b),&
            integrate_pchip(x, im(y), a, b) )
  end function

  function interpolate_pchip(x, y, p) result(r)
    ! This function constructs a piecewise hermitian cubic interpolation of an array y(x) based on discrete numerical data,
    ! and evaluates the interpolation at points p. Note that the mesh spacing of x does not necessarily have to be uniform.
    real(wp), intent(in)  :: x(:)
    real(wp), intent(in)  :: y(:)
    real(wp), intent(in)  :: p(:)
    real(wp), allocatable :: r(:)
    real(wp), allocatable :: d(:)
    integer               :: err
    integer               :: n
    integer               :: m

    ! Check the size of input arrays
    n = min(size(x), size(y))
    m = size(p)

    ! Allocate memory for the derivatives and results
    allocate(d(n))
    allocate(r(m))

    ! Create a PCHIP interpolation of the input data
    call dpchez(n, x, y, d, .false., 0, 0, err)

    ! Extract the interpolated data at provided points
    call dpchfe(n, x, y, d, 1, .false., m, p, r, err)

    ! Clean up workspace memory
    deallocate(d)
  end function

  function interpolate_pchip_cx(x, y, p) result(r)
    ! Wrapper for interpolate_pchip that accepts complex arguments.
    real(wp),    intent(in)  :: x(:)
    real(wp),    intent(in)  :: p(:)
    complex(wp), intent(in)  :: y(:)
    complex(wp), allocatable :: r(:)

    ! Interpolate the real and imaginary parts separately
    r = cx( interpolate_pchip(x, re(y), p),&
            interpolate_pchip(x, im(y), p) )
  end function

  !---------------------------------------------------------------------------!
  !                          MATH PARSING PROCEDURES                          !
  !---------------------------------------------------------------------------!

  impure function evaluate_scalar(expression, location) result(field)
    !! This function takes a mathematical function of z and a position array
    !! as input, and returns the function evaluated at each of the locations
    !! in the array. This is useful for initializing a numeric array of real
    !! values from a corresponding analytical expression.
    use fparser

    character(*), intent(in) :: expression
    real(wp),     intent(in) :: location(:)
    real(wp),    allocatable :: field(:)
    integer                  :: n

    ! Make sure the expression is non-empty
    if (scan(expression, '0123456789piz') <= 0) then
      call error('Invalid scalar expression: "' // expression // '"')
    end if
    
    ! Initialize the function parser
    call initf(1)
    call parsef(1, expression, ['pi', 'z '])

    ! Allocate memory for the output
    allocate(field(size(location)))

    ! Evaluate the parsed function
    do n=1,size(location)
      field(n) = evalf(1, [pi, location(n)])
    end do
  end function

  impure function evaluate_vector(expression, location) result(field)
    !! This function takes a mathematical function of z and a position array
    !! as input, and returns the function evaluated at each of the locations
    !! in the array. This is useful for initializing a numeric array of real
    !! vectors from a corresponding analytical expression.
    use fparser

    character(*), intent(in) :: expression
    real(wp),     intent(in) :: location(:)
    real(wp),    allocatable :: field(:,:)
    integer                  :: n, sep(4)

    ! Allocate memory for the output
    allocate(field(3,size(location)))

    ! Find the vector delimiters
    sep(1) = scan(expression, '[', back=.false.)
    sep(2) = scan(expression, ',', back=.false.)
    sep(3) = scan(expression, ',', back=.true. )
    sep(4) = scan(expression, ']', back=.true. )

    ! Make sure the expressions are non-empty
    if (sep(1) <= 0 .or. any(sep(2:4)-sep(1:3) <= 1)) then
      call error('Invalid vector expression: "' // expression // '"')
    end if
    if (scan(expression(sep(1)+1:sep(2)-1), '0123456789piz') <= 0 .or. &
        scan(expression(sep(2)+1:sep(3)-1), '0123456789piz') <= 0 .or. &
        scan(expression(sep(3)+1:sep(4)-1), '0123456789piz') <= 0) then
      call error('Invalid vector expression: "' // expression // '"')
    end if
    
    ! Initialize the function parser
    call initf(3)
    call parsef(1, expression(sep(1)+1:sep(2)-1), ['pi', 'z '])
    call parsef(2, expression(sep(2)+1:sep(3)-1), ['pi', 'z '])
    call parsef(3, expression(sep(3)+1:sep(4)-1), ['pi', 'z '])

    ! Evaluate the parsed function
    do n=1,size(location)
      field(1,n) = evalf(1, [pi, location(n)])
      field(2,n) = evalf(2, [pi, location(n)])
      field(3,n) = evalf(3, [pi, location(n)])
    end do
  end function
end module
