!> Author:   Jabir Ali Ouassou
!> Category: Foundation
!>
!> This module implements standard operations from linear algebra, including
!> including matrix construction, inversion, commutation, and taking traces.

module matrix_m
    use :: basic_m
contains
    pure function identity(n) result(R)
    !!  Constructs an n×n identity matrix.
        integer,  intent(in)      :: n !! Matrix dimension
        real(wp), dimension(n, n) :: R !! Identity matrix [n×n]

        integer :: i, j

        ! Initialize by exploiting integer arithmetic to avoid multiple passes
        do i = 1, n
            do j = 1, n
                R(j, i) = (i/j)*(j/i)
            end do
        end do
    end function

    pure function matrix_inverse_re(A) result(R)
    !!  Wrapper for matrix_inverse_cx that operates on real matrices.
        real(wp), dimension(:, :), intent(in)       :: A !! Matrix A [n×n]
        real(wp), dimension(size(A, 1), size(A, 1)) :: R !! Matrix R=A¯¹

        R = re(matrix_inverse_cx(cx(A)))
    end function

    pure function matrix_inverse_cx(A) result(R)
    !!  Inverts a square matrix via Gauss-Jordan elimination with partial pivot
    !!  (general) or a cofactoring algorithm (2x2 matrices). The implementation
    !!  is based on Algorithm #2 in "Efficient matrix inversion via Gauss-Jordan
    !!  elimination and its parallelization" by E.S. Quintana et al. (1998).
        complex(wp), dimension(:, :), intent(in)       :: A !! Matrix A [n×n]
        complex(wp), dimension(size(A, 1), size(A, 1)) :: R !! Matrix R=A¯¹
        integer,     dimension(size(A, 1))             :: P

        complex(wp) :: Q
        integer     :: i, j

        select case (size(A, 1))
            case (1)
                ! Trivial case
                R(1, 1) = 1/A(1, 1)

            case (2)
                ! Inverse determinant
                Q = 1/(A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1))

                ! Inverse matrix
                R(1, 1) = +Q*A(2, 2)
                R(2, 1) = -Q*A(2, 1)
                R(1, 2) = -Q*A(1, 2)
                R(2, 2) = +Q*A(1, 1)

            case default
                ! Permutation array
                P = [(i, i=1, size(A, 1))]

                ! Matrix copy
                R = A

                ! Matrix inversion
                do i = 1, size(A, 1)
                    ! Pivoting procedure
                    j = (i - 1) + maxloc(abs(A(i:, i)), 1)
                    P([i, j])    = P([j, i])
                    R([i, j], :) = R([j, i], :)

                    ! Jordan transformation
                    Q = R(i, i)
                    R(:, i) = [R(:i-1, i), (0.0_wp, 0.0_wp), R(i+1:, i)]/(-Q)

                    R = R + matmul(R(:, [i]), R([i], :))
                    R(i, :) = [R(i, :i-1), (1.0_wp, 0.0_wp), R(i, i+1:)]/(+Q)
                end do

                ! Pivot inversion
                R(:, P) = R
        end select
    end function

    pure function matrix_trace(A) result(r)
    !!  Calculates the trace of a general complex matrix.
        complex(wp), dimension(:, :), intent(in) :: A !! Matrix [n×m]
        complex(wp)                              :: r !! r = Tr(A)

        integer :: n

        r = 0
        do n = 1, min(size(A, 1), size(A, 2))
            r = r + A(n, n)
        end do
    end function

    pure function commutator(A, B) result(R)
    !!  Calculates the commutator between two complex square matrices.
        complex(wp), dimension(:, :),                   intent(in) :: A !! Left  matrix [n×n]
        complex(wp), dimension(size(A, 1), size(A, 1)), intent(in) :: B !! Right matrix [n×n]
        complex(wp), dimension(size(A, 1), size(A, 1))             :: R !! Commutator R = [A,B]

        R = matmul(A, B) - matmul(B, A)
    end function

    pure function anticommutator(A, B) result(R)
    !!  Calculates the anticommutator between two complex square matrices.
        complex(wp), dimension(:, :),                   intent(in) :: A !! Left  matrix [n×n]
        complex(wp), dimension(size(A, 1), size(A, 1)), intent(in) :: B !! Right matrix [n×n]
        complex(wp), dimension(size(A, 1), size(A, 1))             :: R !! Anticommutator R = {A,B}

        R = matmul(A, B) + matmul(B, A)
    end function

    pure function vector_diag(A) result(r)
    !!  Extracts the diagonal of a general complex matrix.
        complex(wp), dimension(:, :), intent(in)            :: A !! Matrix [n×m]
        complex(wp), dimension(min(size(A, 1), size(A, 2))) :: r !! r = Diag(A)

        integer :: n

        do n = 1, size(r)
            r(n) = A(n, n)
        end do
    end function

    pure function matrix_diag(A, B) result(R)
    !!  Constructs a block-diagonal matrix R from two general matrices A and B.
        complex(wp), dimension(:, :), intent(in) :: A !! Left  matrix [n×m]
        complex(wp), dimension(:, :), intent(in) :: B !! Right matrix [p×q]
        complex(wp), dimension(size(A, 1) + size(B, 1), size(A, 2) + size(B, 2)) :: R !! R = Diag(A,B)

        R = 0.0_wp
        R(:size(A, 1), :size(A, 2)) = A
        R(size(A, 1) + 1:, size(A, 2) + 1:) = B
    end function
end module
