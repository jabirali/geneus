!> Author:   Jabir Ali Ouassou
!> Date:     2015-04-26
!> Category: Foundation
!>
!> This file defines functions that perform some common matrix operations.

module matrix_m
  use :: math_m
  private

  ! Declare which routines to export
  public :: trace, diag, inverse, identity, commutator, anticommutator

  ! Define common identity matrices
  real(wp), parameter, public :: identity2(2,2) = reshape([1,0,0,1],[2,2])                         !! 2×2 identity matrix
  real(wp), parameter, public :: identity3(3,3) = reshape([1,0,0,0,1,0,0,0,1],[3,3])               !! 3×3 identity matrix
  real(wp), parameter, public :: identity4(4,4) = reshape([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],[4,4]) !! 4×4 identity matrix

  ! Declare public interfaces
  interface trace
    !! Public interface for functions that calculate a matrix trace.
    module procedure matrix_trace
  end interface

  interface inverse
    !! Public interface for functions that calculate a matrix inverse.
    module procedure matrix_inverse
  end interface
  
  interface diag
    !! Public interface for functions that either extract a matrix diagonal, or construct a diagonal matrix
    module procedure matrix_diag, vector_diag
  end interface
contains
  pure function identity(n) result(A)
    !! Constructs an N×N identity matrix.
    integer, intent(in) :: n        !! Matrix dimension
    real(wp)            :: A(n,n)   !! Identity matrix
    integer             :: i, j

    ! Initialize by exploiting integer arithmetic to avoid multiple passes
    do i=1,n
      do j=1,n
        A(j,i) = (i/j)*(j/i)
      end do
    end do
  end function

  pure function matrix_inverse(A) result(B)
    !! Performs a direct calculation of the inverse of an N×N matrix, where N≤4.
    complex(wp), intent(in) :: A(:,:)                   !! Matrix
    complex(wp)             :: B(size(A,1),size(A,2))   !! Inverse matrix

    select case(size(A))
      case(2**2)
        B = matrix_inverse2(A)
      case(3**2)
        B = matrix_inverse3(A)
      case(4**2)
        B = matrix_inverse4(A)
      case default
        B = 0.0_wp
    end select
  end function

  pure function matrix_inverse2(A) result(B)
    !! Performs a direct calculation of the inverse of a 2×2 matrix.
    complex(wp), intent(in) :: A(2,2)   !! Matrix
    complex(wp)             :: B(2,2)   !! Inverse matrix
    complex(wp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)
  end function

  pure function matrix_inverse3(A) result(B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    !! [Based on the subroutine M33INV by David G. Simpson, NASA.]
    complex(wp), intent(in) :: A(3,3)   !! Matrix
    complex(wp)             :: B(3,3)   !! Inverse matrix
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

  pure function matrix_inverse4(A) result(B)
    !! Performs a direct calculation of the inverse of a 4×4 matrix.
    !! [Based on the subroutine M44INV by David G. Simpson, NASA.]
    complex(wp), intent(in) :: A(4,4)   !! Matrix
    complex(wp)             :: B(4,4)   !! Inverse matrix
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

  pure function matrix_trace(A) result(r)
    !! Calculate the trace of a complex matrix.
    complex(wp), intent(in)  :: A(:,:)   !! Matrix
    complex(wp)              :: r        !! Tr(A)
    integer                  :: n

    r = 0
    do n = 1,min(size(A,1),size(A,2))
      r = r + A(n,n)
    end do
  end function

  pure subroutine factorize(a,p)
    !! Perform an in-place PLU decomposition of a square matrix A, similar to LAPACK's *gebtrf subroutines.
    !! [Based on the example code presented on http://rosettacode.org/wiki/LU_decomposition#Fortran]
    !! @TODO: This function is still experimental, and needs testing before being used for serious stuff.
    complex(wp), intent(inout) :: a(:,:)
    integer,     intent(  out) :: p(size(a,1))
    integer                    :: n, m, i, j

    ! Check the size of the input matrix
    n = size(a,1)

    ! Initialize the permutation vector
    p = [ ( i, i=1,n ) ]

    ! Perform the PLU factorization
    do i = 1,n-1
      ! Permutation
      m = maxloc(abs(a(p(i:),i)),1) + i-1
      if (m /= i) then
        p([i, m]) = p([m, i])
      end if

      ! Factorizatrion
      a(p(i+1:),i) = a(p(i+1:),i) / a(p(i),i) 
      do j=i+1,n
        a(p(i+1:),j) = a(p(i+1:),j) - a(p(i+1:),i) * a(p(i),j)
      end do
    end do
  end subroutine

  pure function commutator(A, B) result(C)
    !! Calculate the commutator between two complex square matrices of the same dimension.
    complex(wp), intent(in) :: A(:,:)                   !! Left  matrix
    complex(wp), intent(in) :: B(size(A,1),size(A,1))   !! Right matrix
    complex(wp)             :: C(size(A,1),size(A,1))   !! Commutator [A,B]

    C = matmul(A,B) - matmul(B,A)
  end function

  pure function anticommutator(A, B) result(C)
    !! Calculate the anticommutator between two complex square matrices of the same dimension.
    complex(wp), intent(in)  :: A(:,:)                   !! Left  matrix
    complex(wp), intent(in)  :: B(size(A,1),size(A,1))   !! Right matrix
    complex(wp)              :: C(size(A,1),size(A,1))   !! Anticommutator {A,B}

    C = matmul(A,B) + matmul(B,A)
  end function

  pure function vector_diag(A) result(r)
    !! Extract the diagonal of a complex matrix.
    complex(wp), intent(in)  :: A(:,:)                        !! Matrix
    complex(wp)              :: r(min(size(A,1),size(A,2)))   !! Diag(A)
    integer                  :: n

    do n = 1,size(r)
      r(n) = A(n,n)
    end do
  end function

  pure function matrix_diag(A,B) result(C)
    !! Construct a block-diagonal matrix C from two matrices A and B.
    complex(wp), intent(in) :: A(:,:), B(:,:)
    complex(wp)             :: C(size(A,1)+size(B,1), size(A,2)+size(B,2))

    C = 0.0_wp
    C(:size(A,1), :size(A,2))     = A
    C(size(A,1)+1:, size(A,2)+1:) = B
  end function
end module
