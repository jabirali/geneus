! This file defines a module containing the machine size of single-precision, double-precision, and quadruple-precision
! floating point numbers; to declare the floating point precision of a variable, use real(sp), real(dp), or real(qp) as
! the type of the variable. It also defines some common mathematical functions, like matrix inversion of small matrices.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-09-23
! Updated: 2015-09-23

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

contains
  pure function matinv2(A) result(Ainv)
    ! Performs a direct calculation of the inverse of a 2×2 matrix.
    complex(dp), intent(in) :: A(2,2)
    complex(dp)             :: Ainv(2,2)
    complex(dp)             :: detinv

    ! Calculate the inverse determinant of this matrix
    detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse matrix as a result
    Ainv(1,1) = +detinv * A(2,2)
    Ainv(1,2) = -detinv * A(1,2)
    Ainv(2,1) = -detinv * A(2,1)
    Ainv(2,2) = +detinv * A(1,1)
  end function

  pure function matinv3(A) result(Ainv)
    ! Performs a direct calculation of the inverse of a 3×3 matrix.
    ! [Based on the subroutine M33INV by David G. Simpson, NASA.]
    complex(dp), intent(in) :: A(3,3)
    complex(dp)             :: Ainv(3,3)
    complex(dp)             :: cofactor(3,3)
    complex(dp)             :: det

    ! Calculate the determinant of the matrix
    det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) &
        - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
        + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

    ! Calculate the cofactors of the matrix
    cofactor(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
    cofactor(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    cofactor(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    cofactor(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
    cofactor(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
    cofactor(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
    cofactor(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
    cofactor(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
    cofactor(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

    ! Calculate the inverse matrix
    Ainv = transpose(cofactor)/det
  end function

  pure function matinv4(A) result(Ainv)
    ! Performs a direct calculation of the inverse of a 4×4 matrix.
    ! [Based on the subroutine M44INV by David G. Simpson, NASA.]
    complex(dp), intent(in) :: A(4,4)
    complex(dp)             :: Ainv(4,4)
    complex(dp)             :: cofactor(4,4)
    complex(dp)             :: det

    ! Calculate the determinant of the matrix
    det = A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))) &
        - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))) &
        + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))) &
        - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))

    ! Calculate the cofactors of the matrix
    cofactor(1,1) = A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))
    cofactor(1,2) = A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))
    cofactor(1,3) = A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
    cofactor(1,4) = A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
    cofactor(2,1) = A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))
    cofactor(2,2) = A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))
    cofactor(2,3) = A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
    cofactor(2,4) = A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
    cofactor(3,1) = A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2))
    cofactor(3,2) = A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3))
    cofactor(3,3) = A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1))
    cofactor(3,4) = A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2))
    cofactor(4,1) = A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3))
    cofactor(4,2) = A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    cofactor(4,3) = A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2))
    cofactor(4,4) = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

    ! Calculate the inverse matrix
    Ainv = transpose(cofactor)/det
  end function
end module
