! This file defines a module containing the machine size of single-precision, double-precision, and quadruple-precision
! floating point numbers; to declare the floating point precision of a variable, use real(sp), real(dp), or real(qp) as
! the type of the variable. It also defines the working-precision, which will be the default kind for module procedures.
! As for module procedures, this library defines common utility functions for working with complex numbers and matrices.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-09-23
! Updated: 2016-04-06

module math_m
  use :: iso_fortran_env
  use :: stdio_m
  implicit none
  public

  ! Declare floating-point precisions
  integer,     parameter :: sp  = REAL32              !! Single precision
  integer,     parameter :: dp  = REAL64              !! Double precision
  integer,     parameter :: qp  = REAL128             !! Quadruple precision
  integer,     parameter :: wp  = dp                  !! Working precision

  ! Define common mathematical constants
  real(wp),    parameter :: inf = huge(1.0_wp)        !! Numerical infinity
  real(wp),    parameter :: eps = epsilon(1.0_wp)     !! Numerical infinitesimal
  real(wp),    parameter :: pi  = atan(1.0_wp)*4.0_wp !! Circle constant

  ! Define common identity matrices
  real(wp),    parameter :: mateye2(2,2) = reshape([1,0,0,1],[2,2])                         !! 2×2 identity matrix
  real(wp),    parameter :: mateye3(3,3) = reshape([1,0,0,0,1,0,0,0,1],[3,3])               !! 3×3 identity matrix
  real(wp),    parameter :: mateye4(4,4) = reshape([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],[4,4]) !! 4×4 identity matrix

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

  interface evaluate
    module procedure evaluate_scalar_value, evaluate_scalar_field, &
                     evaluate_vector_value, evaluate_vector_field, &
                     evaluate_logical_value
  end interface
contains

  !---------------------------------------------------------------------------!
  !                        COMPLEX NUMBER PROCEDURES                          !
  !---------------------------------------------------------------------------!

  pure elemental function re(z) result(x)
    !! Returns the real part of a complex number z=x+iy.
    complex(wp), intent(in) :: z
    real(wp)                :: x

    x = real(z,kind=wp)
  end function

  pure elemental function im(z) result(y)
    !! Returns the imaginary part of a complex number z=x+iy.
    complex(wp), intent(in) :: z
    real(wp)                :: y

    y = aimag(z)
  end function

  pure elemental function cx(x,y) result(z)
    !! Returns the complex number z=x+iy.
    real(wp),           intent(in) :: x
    real(wp), optional, intent(in) :: y
    complex(wp)                    :: z

    if (present(y)) then
      z = cmplx(x,y,kind=wp)
    else
      z = cmplx(x,0,kind=wp)
    end if
  end function

  !---------------------------------------------------------------------------!
  !                       ELEMENTARY VECTOR PROCEDURES                        !
  !---------------------------------------------------------------------------!

  pure subroutine linspace(array, first, last)
    !! Populates an array with elements from 'first' to 'last', inclusive.
    real(wp), intent(inout) :: array(:)
    real(wp), intent(in   ) :: first, last
    integer                 :: n

    do n=1,size(array)
      array(n) = first + ((last-first)*(n-1))/(size(array)-1)
    end do
  end subroutine

  pure function unitvector(v) result(r)
    !! Rescales a vector to a unit vector. If the input vector has zero length,
    !! then the output vector will be a zero vector instead of a unit vector.
    real(wp), intent(in) :: v(3)
    real(wp)             :: r(3)

    r = v/(norm2(v)+eps)
  end function

  !---------------------------------------------------------------------------!
  !                       ELEMENTARY MATRIX PROCEDURES                        !
  !---------------------------------------------------------------------------!

  pure function mateye(n) result(A)
    !! Constructs an n×n identity matrix.
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

  pure function matinv2(A) result(B)
    !! Performs a direct calculation of the inverse of a 2×2 matrix.
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
    !! Performs a direct calculation of a 2×2 matrix left division.
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
    !! Performs a direct calculation of a 2×2 matrix right division.
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
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    !! [Based on the subroutine M33INV by David G. Simpson, NASA.]
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
    !! Performs a direct calculation of the inverse of a 4×4 matrix.
    !! [Based on the subroutine M44INV by David G. Simpson, NASA.]
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
  !                        ADVANCED MATRIX PROCEDURES                         !
  !---------------------------------------------------------------------------!

  pure function commutator(A, B) result(C)
    !! Calculate the commutator between two complex square matrices of the same dimension.
    complex(wp), intent(in)  :: A(:,:)
    complex(wp), intent(in)  :: B(:,:)
    complex(wp), allocatable :: C(:,:)

    C = matmul(A,B) - matmul(B,A)
  end function

  pure function anticommutator(A, B) result(C)
    !! Calculate the anticommutator between two complex square matrices of the same dimension.
    complex(wp), intent(in)  :: A(:,:)
    complex(wp), intent(in)  :: B(:,:)
    complex(wp), allocatable :: C(:,:)

    C = matmul(A,B) + matmul(B,A)
  end function

  pure function trace(A) result(r)
    !! Calculate the trace of a complex matrix.
    complex(wp), intent(in)  :: A(:,:)
    complex(wp)              :: r
    integer                  :: n

    r = 0
    do n = 1,min(size(A,1),size(A,2))
      r = r + A(n,n)
    end do
  end function

  pure function diag(A) result(r)
    !! Extract the diagonal of a complex matrix.
    complex(wp), intent(in)  :: A(:,:)
    complex(wp), allocatable :: r(:)
    integer                  :: n

    allocate(r(min(size(A,1),size(A,2))))
    do n = 1,size(r)
      r(n) = A(n,n)
    end do
  end function

  !---------------------------------------------------------------------------!
  !                       ELEMENTARY CALCULUS PROCEDURES                      !
  !---------------------------------------------------------------------------!

  pure function differentiate_linear(x, y) result(r)
    !! This function calculates the numerical derivative of an array y with respect to x, using a central difference approximation
    !! at the interior points and forward/backward difference approximations at the exterior points. Note that since all the three
    !! approaches yield two-point approximations of the derivative, the mesh spacing of x does not necessarily have to be uniform.
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
    !! Complex version of differentiate_linear.
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
    !! This function calculates the integral of an array y with respect to x using a trapezoid
    !! approximation. Note that the mesh spacing of x does not necessarily have to be uniform.
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
    !! Complex version of integrate_linear.
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
    !! This function constructs a piecewise hermitian cubic interpolation of an array y(x) based on
    !! discrete numerical data, and subsequently evaluates the integral of the interpolation in the
    !! range (a,b). Note that the mesh spacing of x does not necessarily have to be uniform.
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
    !! Wrapper for integrate_pchip that accepts complex arguments.
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
    !! This function constructs a piecewise hermitian cubic interpolation of an array y(x) based on discrete numerical data,
    !! and evaluates the interpolation at points p. Note that the mesh spacing of x does not necessarily have to be uniform.
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
    !! Wrapper for interpolate_pchip that accepts complex arguments.
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

  impure subroutine evaluate_logical_value(expression, value)
    !! This subroutine takes a scalar logical expression as input, and returns the value.
    !!
    !! Usage:
    !! 
    !!     call evaluate_logical_value('F', output)
    !!     call evaluate_logical_value('T', output)
    !!

    character(*), intent(in)  :: expression !! Either the character 'T' or 'F'
    logical,      intent(out) :: value      !! Result of parsing the expression

    select case(expression)
      case ('T', 't')
        value = .true.
      case ('F', 'f')
        value = .false.
      case default
        call error('Invalid logical expression: "' // trim(expression) // '"')
    end select
  end subroutine

  impure subroutine evaluate_scalar_value(expression, value)
    !! This subroutine takes a scalar mathematical expression as input, and returns the value.
    !!
    !! Usage:
    !! 
    !!     call evaluate_scalar_value('0',                    output)
    !!     call evaluate_scalar_value('sin(0.3*pi)*exp(-pi)', output)
    !!
    use :: fparser

    character(*), intent(in)  :: expression  !! Scalar-valued mathematical expression
    real(wp),     intent(out) :: value       !! Result of parsing the expression

    ! Make sure the expression is non-empty
    if (scan(expression, '0123456789pi') <= 0) then
      call error('Invalid scalar expression: "' // trim(expression) // '"')
    end if

    ! Initialize the function parser
    call initf(1)
    call parsef(1, expression, ['pi'])

    ! Evaluate the parsed function
    value = evalf(1, [pi])
  end subroutine

  impure subroutine evaluate_scalar_field(expression, domain, value)
    !! This subroutine takes a scalar mathematical function of some variable 'z' as input,  along 
    !! with an array with discrete values for that variable 'z'.  It parses the provided function,
    !! evaluates it at each 'z'-value in the array, and then returns the discretized scalar field.
    !!
    !! Usage:
    !! 
    !!     call evaluate_scalar_value('0',                    input(1:n), output(1:n))
    !!     call evaluate_scalar_value('sin(pi*z)*exp(-pi*z)', input(1:n), output(1:n))
    !!
    use :: fparser

    character(*), intent(in) :: expression   !! Scalar-valued function of position 'z'
    real(wp),     intent(in) :: domain(:)    !! Domain of the independent variable 'z'
    real(wp),    allocatable :: value(:)     !! Result of evaluating the field at each point of the domain
    integer                  :: n

    ! Make sure the expression is non-empty
    if (scan(expression, '0123456789piz') <= 0) then
      call error('Invalid scalar expression: "' // trim(expression) // '"')
    end if

    ! Initialize the function parser
    call initf(1)
    call parsef(1, expression, ['pi', 'z '])

    ! Allocate memory for the output
    allocate(value(size(domain)))

    ! Evaluate the parsed function
    do n=1,size(domain)
      value(n) = evalf(1, [pi, domain(n)])
    end do
  end subroutine

  impure subroutine evaluate_vector_value(expression, value)
    !! This subroutine takes a vector mathematical expression as input, and returns the value.
    !!
    !! Usage:
    !! 
    !!     call evaluate_scalar_value('[0,0,0]',                     output(1:3))
    !!     call evaluate_scalar_value('[sin(0.3*pi),0,cos(0,3*pi)]', output(1:3))
    !!
    use :: fparser

    character(*), intent(in)  :: expression  !! Vector-valued mathematical expression
    real(wp),     intent(out) :: value(3)    !! Result of parsing the expression
    integer                   :: sep(4)

    ! Find the vector delimiters
    sep(1) = scan(expression, '[', back=.false.)
    sep(2) = scan(expression, ',', back=.false.)
    sep(3) = scan(expression, ',', back=.true. )
    sep(4) = scan(expression, ']', back=.true. )

    ! Make sure the expressions are non-empty
    if (sep(1) <= 0 .or. any(sep(2:4)-sep(1:3) <= 1)) then
      call error('Invalid vector expression: "' // trim(expression) // '"')
    end if
    if (scan(expression(sep(1)+1:sep(2)-1), '0123456789pi') <= 0 .or. &
        scan(expression(sep(2)+1:sep(3)-1), '0123456789pi') <= 0 .or. &
        scan(expression(sep(3)+1:sep(4)-1), '0123456789pi') <= 0) then
      call error('Invalid vector expression: "' // trim(expression) // '"')
    end if

    ! Initialize the function parser
    call initf(3)
    call parsef(1, expression(sep(1)+1:sep(2)-1), ['pi'])
    call parsef(2, expression(sep(2)+1:sep(3)-1), ['pi'])
    call parsef(3, expression(sep(3)+1:sep(4)-1), ['pi'])

    ! Evaluate the parsed function
    value(1) = evalf(1, [pi])
    value(2) = evalf(2, [pi])
    value(3) = evalf(3, [pi])
  end subroutine

  impure subroutine evaluate_vector_field(expression, domain, value)
    !! This subroutine takes a vector mathematical function of some variable 'z' as input,  along 
    !! with an array with discrete values for that variable 'z'.  It parses the provided function,
    !! evaluates it at each 'z'-value in the array, and then returns the discretized scalar field.
    !!
    !! Usage:
    !! 
    !!     call evaluate_scalar_value('[0,0,0]',                     input(1:n), output(1:3,1:n))
    !!     call evaluate_scalar_value('[sin(pi*z/2),0,cos(pi*z/2)]', input(1:n), output(1:3,1:n))
    !!
    use :: fparser

    character(*), intent(in) :: expression   !! Vector-valued function of position 'z'
    real(wp),     intent(in) :: domain(:)    !! Domain of the independent variable 'z'
    real(wp),    allocatable :: value(:,:)   !! Result of evaluating the field at each point of the domain
    integer                  :: n, sep(4)

    ! Allocate memory for the output
    allocate(value(3,size(domain)))

    ! Find the vector delimiters
    sep(1) = scan(expression, '[', back=.false.)
    sep(2) = scan(expression, ',', back=.false.)
    sep(3) = scan(expression, ',', back=.true. )
    sep(4) = scan(expression, ']', back=.true. )

    ! Make sure the expressions are non-empty
    if (sep(1) <= 0 .or. any(sep(2:4)-sep(1:3) <= 1)) then
      call error('Invalid vector expression: "' // trim(expression) // '"')
    end if
    if (scan(expression(sep(1)+1:sep(2)-1), '0123456789piz') <= 0 .or. &
        scan(expression(sep(2)+1:sep(3)-1), '0123456789piz') <= 0 .or. &
        scan(expression(sep(3)+1:sep(4)-1), '0123456789piz') <= 0) then
      call error('Invalid vector expression: "' // trim(expression) // '"')
    end if

    ! Initialize the function parser
    call initf(3)
    call parsef(1, expression(sep(1)+1:sep(2)-1), ['pi', 'z '])
    call parsef(2, expression(sep(2)+1:sep(3)-1), ['pi', 'z '])
    call parsef(3, expression(sep(3)+1:sep(4)-1), ['pi', 'z '])

    ! Evaluate the parsed function
    do n=1,size(domain)
      value(1,n) = evalf(1, [pi, domain(n)])
      value(2,n) = evalf(2, [pi, domain(n)])
      value(3,n) = evalf(3, [pi, domain(n)])
    end do
  end subroutine
end module
