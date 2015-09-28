! This file defines a module containing the machine size of single-precision, double-precision, and quadruple-precision
! floating point numbers; to declare the floating point precision of a variable, use real(sp), real(dp), or real(qp) as
! the type of the variable. It also defines some common mathematical functions, like matrix inversion of small matrices.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-09-23
! Updated: 2015-09-28

module mod_math
  use, intrinsic :: iso_fortran_env

  ! Declare floating-point precisions
  integer,     parameter :: sp  = REAL32
  integer,     parameter :: dp  = REAL64
  integer,     parameter :: qp  = REAL128

  ! Define common mathematical constants
  real(dp),    parameter :: inf = huge(1.0_dp)
  real(dp),    parameter :: pi  = atan(1.0_dp)*4.0_dp
  complex(dp), parameter :: i   = (0.0_dp,1.0_dp)

  ! Define common identity matrices
  real(dp),    parameter :: mateye2(2,2) = reshape([1,0,0,1],[2,2])
  real(dp),    parameter :: mateye3(3,3) = reshape([1,0,0,0,1,0,0,0,1],[3,3])
  real(dp),    parameter :: mateye4(4,4) = reshape([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],[4,4])

contains
  pure function mateye(n) result(A)
    ! Constructs an n×n identity matrix.
    integer, intent(in) :: n
    integer             :: i, j
    real(dp)            :: A(n,n)

    ! Initialize by exploiting integer arithmetic to avoid multiple passes
    do i=1,n
      do j=1,n
        A(i,j) = (i/j)*(j/i)
      end do
    end do
  end function

  function matrnd(n) result(A)
    ! Constructs an n×n random matrix.
    integer, allocatable :: seed(:)
    integer              :: n, m, u
    complex(dp)          :: A(n,n)
    real(dp)             :: R(n,n)
    real(dp)             :: I(n,n)

    ! Find the minimum size of an RNG seed
    call random_seed(size = m)

    ! Allocate memory for the RNG seed
    allocate(seed(m))

    ! Initialize the seed using random numbers from the operating system
    open(newunit=u, file="/dev/urandom", access="stream", form="unformatted", action="read", status="old")
    read(unit=u) seed
    close(unit=u)

    ! Initialize the RNG using the seed
    call random_seed(put = seed)

    ! Construct a random complex matrix
    call random_number(R)
    call random_number(I)
    A = cmplx(R,I,kind=dp)

    ! Deallocate dynamic memory
    deallocate(seed)
  end function

  pure function matinv2(A) result(B)
    ! Performs a direct calculation of the inverse of a 2×2 matrix.
    complex(dp), intent(in) :: A(2,2)
    complex(dp)             :: B(2,2)
    complex(dp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)
  end function

  pure function matinv3(A) result(B)
    ! Performs a direct calculation of the inverse of a 3×3 matrix.
    ! [Based on the subroutine M33INV by David G. Simpson, NASA.]
    complex(dp), intent(in) :: A(3,3)
    complex(dp)             :: B(3,3)
    complex(dp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3)-A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3)-A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2)-A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3)-A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3)-A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2)-A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3)-A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3)-A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2)-A(1,2)*A(2,1))
  end function

  pure function matinv4(A) result(B)
    ! Performs a direct calculation of the inverse of a 4×4 matrix.
    ! [Based on the subroutine M44INV by David G. Simpson, NASA.]
    complex(dp), intent(in) :: A(4,4)
    complex(dp)             :: B(4,4)
    complex(dp)             :: detinv

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

  pure function matdivl2(A,B) result(C)
    ! Performs a direct calculation of a 2×2 matrix left division.
    complex(dp), intent(in) :: A(2,2)
    complex(dp), intent(in) :: B(2,2)
    complex(dp)             :: C(2,2)
    complex(dp)             :: detinv

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
    complex(dp), intent(in) :: A(2,2)
    complex(dp), intent(in) :: B(2,2)
    complex(dp)             :: C(2,2)
    complex(dp)             :: detinv

    ! Calculate the inverse determinant of the right matrix
    detinv = 1/(B(1,1)*B(2,2) - B(1,2)*B(2,1))

    ! Calculate the elements of the resulting matrix
    C(1,1) = detinv * (A(1,1)*B(2,2) - A(1,2)*B(2,1))
    C(2,1) = detinv * (A(2,1)*B(2,2) - A(2,2)*B(2,1))
    C(1,2) = detinv * (A(1,2)*B(1,1) - A(1,1)*B(1,2))
    C(2,2) = detinv * (A(2,2)*B(1,1) - A(2,1)*B(1,2))
  end function
end module
