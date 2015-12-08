! This file defines a module containing the machine size of single-precision, double-precision, and quadruple-precision
! floating point numbers; to declare the floating point precision of a variable, use real(sp), real(dp), or real(qp) as
! the type of the variable. It also defines the working-precision, which will be the default kind for module procedures.
! As for module procedures, this library defines common utility functions for working with complex numbers and matrices.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-09-23
! Updated: 2015-10-04

module mod_math
  use :: iso_fortran_env
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
contains

  !----------------------------------------------------------------------------!
  !                        COMPLEX NUMBER PROCEDURES                           !
  !----------------------------------------------------------------------------!

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

  !----------------------------------------------------------------------------!
  !                       ELEMENTARY MATRIX PROCEDURES                         !
  !----------------------------------------------------------------------------!

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

  !----------------------------------------------------------------------------!
  !                       ELEMENTARY CALCULUS PROCEDURES                       !
  !----------------------------------------------------------------------------!

  pure function differentiate(x, y) result(r)
    ! This function calculates the numerical derivative of an array y with respect to x, using a central difference approximation 
    ! at the interior points and forward/backward difference approximations at the exterior points. Note that since all the three
    ! approaches yield two-point approximations of the derivative, the mesh spacing of x does not necessarily have to be uniform.
    real(wp), intent(in)  :: x(:)
    real(wp), intent(in)  :: y(:)
    real(wp), allocatable :: r(:)
    integer               :: n

    ! Check the size of input arrays
    n = max(size(x), size(y))

    ! Allocate memory for the results
    allocate(r(n))

    ! Differentiate using finite differences
    r(   1   ) = (y( 1+1 ) - y(  1  ))/(x( 1+1 ) - x(  1  ))
    r(1+1:n-1) = (y(1+2:n) - y(1:n-2))/(x(1+2:n) - x(1:n-2))
    r(   n   ) = (y(  n  ) - y( n-1 ))/(x(  n  ) - x( n-1 ))
  end function

  pure function integrate(x, y) result(r)
    ! This function calculates the integral of an array y with respect to x using a trapezoid 
    ! approximation. Note that the mesh spacing of x does not necessarily have to be uniform.
    real(wp), intent(in)  :: x(:)
    real(wp), intent(in)  :: y(:)
    real(wp)              :: r
    integer               :: n

    ! Check the size of input arrays
    n = max(size(x), size(y))

    ! Integrate using the trapezoidal rule
    r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
  end function

  !pure function interpolate(x, y, p) result(r)
  !end function

end module
