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
    !! Invert a general N×N matrix using Gauss-Jordan elimination with partial pivoting.
    !! In the special case N=2, the inverse is evaluated using a cofactoring algorithm.
    !! [This implementation is based on Algorithm #2 in "Efficient matrix inversion via 
    !! Gauss-Jordan elimination and its parallelization" by E.S. Quintana et al. (1998)]
    complex(wp), dimension(:,:), intent(in)     :: A
    complex(wp), dimension(size(A,1),size(A,1)) :: B
    integer,     dimension(size(A,1))           :: P
    complex(wp)                                 :: Q
    integer                                     :: i, j

    select case (size(A,1))
      case (1)
        ! Trivial case
        B(1,1) = 1/A(1,1)

      case (2)
        ! Inverse determinant
        Q = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

        ! Inverse matrix
        B(1,1) = +Q * A(2,2)
        B(2,1) = -Q * A(2,1)
        B(1,2) = -Q * A(1,2)
        B(2,2) = +Q * A(1,1)

      case default
        ! Permutation array
        P = [ ( i, i=1,size(A,1) ) ]

        ! Matrix copy
        B = A

        ! Matrix inversion
        do i=1,size(A,1)
          ! Pivoting procedure
          j = (i-1) + maxloc(abs(a(i:,i)),1)
          P([i,j])   = P([j,i])
          B([i,j],:) = B([j,i],:)

          ! Jordan transformation
          Q      = B(i,i)
          B(:,i) = [B(:i-1,i), (0.0_wp,0.0_wp), B(i+1:,i)] / (-Q)
          B      = B + matmul(B(:,[i]), B([i],:))
          B(i,:) = [B(i,:i-1), (1.0_wp,0.0_wp), B(i,i+1:)] / (+Q)
        end do

        ! Pivot inversion
        B(:,P) = B
    end select
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
