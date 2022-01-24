!> Author:   Jabir Ali Ouassou
!> Category: Foundation

!> This module defines the type 'spin', representing 2×2 complex matrices in
!> spin space. The module overloads common arithmetic operators to work with
!> the new data type, and defines and exports the Pauli matrices as constants.
!> To make it easier to interact with common differential equation solvers,
!> which often operate on real state vectors, the assignment operator is
!> overloaded to make 'spin' easily importable/exportable to real vectors.

module spin_m
    use :: math_m
    private

    ! Public interface
    public spin
    public inverse, trace, conjg, norm2, sum
    public pauli, pauli0, pauli1, pauli2, pauli3

    ! Type declaration
    type spin
        complex(wp) :: matrix(2, 2) = 0.0_wp
    contains
        ! Overload constructors and operators
        generic :: spin => &
            cons_rscalar, cons_cscalar, &
            cons_cmatrix, cons_rvector, &
            cons_spin
        generic :: assignment(=) => &
            assr_rscalar, assr_cscalar, &
            assr_cmatrix, assr_rvector, &
            assl_cmatrix, assl_rvector
        generic :: operator(+) => &
            addl_rscalar, addr_rscalar, &
            addl_cscalar, addr_cscalar, &
            addl_cmatrix, addr_cmatrix, &
            add_spin
        generic :: operator(-) => &
            subl_rscalar, subr_rscalar, &
            subl_cscalar, subr_cscalar, &
            subl_cmatrix, subr_cmatrix, &
            sub_spin
        generic :: operator(*) => &
            mull_rscalar, mulr_rscalar, &
            mull_cscalar, mulr_cscalar, &
            mull_cmatrix, mulr_cmatrix, &
            mul_spin
        generic :: operator(/) => &
            divr_rscalar, divr_cscalar
        generic :: operator(**) => &
            expr_iscalar

        ! Specific methods for construction
        procedure, nopass, private :: cons_spin    => spin_cons_spin
        procedure, nopass, private :: cons_rscalar => spin_cons_rscalar
        procedure, nopass, private :: cons_cscalar => spin_cons_cscalar
        procedure, nopass, private :: cons_cmatrix => spin_cons_cmatrix
        procedure, nopass, private :: cons_rvector => spin_cons_rvector

        ! Specific implementations of assignments
        procedure, pass(this), private :: assr_rscalar => spin_assr_rscalar
        procedure, pass(this), private :: assr_cscalar => spin_assr_cscalar
        procedure, pass(this), private :: assl_cmatrix => spin_assl_cmatrix
        procedure, pass(this), private :: assr_cmatrix => spin_assr_cmatrix
        procedure, pass(this), private :: assl_rvector => spin_assl_rvector
        procedure, pass(this), private :: assr_rvector => spin_assr_rvector

        ! Specific implementations of addition
        procedure, pass(this), private :: add_spin     => spin_add_spin
        procedure, pass(this), private :: addl_rscalar => spin_addl_rscalar
        procedure, pass(this), private :: addr_rscalar => spin_addr_rscalar
        procedure, pass(this), private :: addl_cscalar => spin_addl_cscalar
        procedure, pass(this), private :: addr_cscalar => spin_addr_cscalar
        procedure, pass(this), private :: addl_cmatrix => spin_addl_cmatrix
        procedure, pass(this), private :: addr_cmatrix => spin_addr_cmatrix

        ! Specific implementations of subtraction
        procedure, pass(this), private :: sub_spin     => spin_sub_spin
        procedure, pass(this), private :: subl_rscalar => spin_subl_rscalar
        procedure, pass(this), private :: subr_rscalar => spin_subr_rscalar
        procedure, pass(this), private :: subl_cscalar => spin_subl_cscalar
        procedure, pass(this), private :: subr_cscalar => spin_subr_cscalar
        procedure, pass(this), private :: subl_cmatrix => spin_subl_cmatrix
        procedure, pass(this), private :: subr_cmatrix => spin_subr_cmatrix

        ! Specific implementations of multiplication
        procedure, pass(this), private :: mul_spin     => spin_mul_spin
        procedure, pass(this), private :: mull_rscalar => spin_mull_rscalar
        procedure, pass(this), private :: mulr_rscalar => spin_mulr_rscalar
        procedure, pass(this), private :: mull_cscalar => spin_mull_cscalar
        procedure, pass(this), private :: mulr_cscalar => spin_mulr_cscalar
        procedure, pass(this), private :: mull_cmatrix => spin_mull_cmatrix
        procedure, pass(this), private :: mulr_cmatrix => spin_mulr_cmatrix

        ! Specific implementations of division
        procedure, pass(this), private :: divr_rscalar => spin_divr_rscalar
        procedure, pass(this), private :: divr_cscalar => spin_divr_cscalar

        ! Specific implementations of exponentiation
        procedure, pass(this), private :: expr_iscalar => spin_expr_iscalar
    end type

    ! Public interfaces
    interface inverse
        module procedure spin_inv
    end interface

    interface trace
        module procedure spin_trace
    end interface

    interface sum
        module procedure spin_sum
    end interface

    interface conjg
        module procedure spin_conjg
    end interface

    interface norm2
        module procedure spin_norm
    end interface

    ! Exported constants
    type(spin), parameter :: pauli0 = &
        spin(reshape([(1, 0), (0, 0),  (0, 0),  (1, 0)], [2, 2], order=[2, 1]))

    type(spin), parameter :: pauli1 = &
        spin(reshape([(0, 0), (1, 0),  (1, 0),  (0, 0)], [2, 2], order=[2, 1]))

    type(spin), parameter :: pauli2 = &
        spin(reshape([(0, 0), (0, -1), (0, 1),  (0, 0)], [2, 2], order=[2, 1]))

    type(spin), parameter :: pauli3 = &
        spin(reshape([(1, 0), (0, 0),  (0, 0), (-1, 0)], [2, 2], order=[2, 1]))

    type(spin), parameter, dimension(0:3) :: pauli = &
        [pauli0, pauli1, pauli2, pauli3]
contains

    !--------------------------------------------------------------------------!
    !                            SPECIFIC CONSTRUCTORS                         !
    !--------------------------------------------------------------------------!

    pure function spin_cons_rscalar(other) result(this)
    !!  Constructs a spin object from a real scalar.
        real(wp), intent(in) :: other
        type(spin)           :: this

        this = other
    end function

    pure function spin_cons_cscalar(other) result(this)
    !!  Constructs a spin object from a complex scalar.
        complex(wp), intent(in) :: other
        type(spin)              :: this

        this = other
    end function

    pure function spin_cons_cmatrix(other) result(this)
    !!  Constructs a spin object from a complex matrix.
        complex(wp), intent(in) :: other(2, 2)
        type(spin)              :: this

        this = other
    end function

    pure function spin_cons_rvector(other) result(this)
    !!  Constructs a spin object from a real vector.
        real(wp), intent(in) :: other(8)
        type(spin)           :: this

        this = other
    end function

    pure function spin_cons_spin(other) result(this)
    !!  Constructs a spin object from an existing one.
        type(spin), intent(in) :: other
        type(spin)             :: this

        this = other
    end function

    !--------------------------------------------------------------------------!
    !                         SPECIFIC IMPORT PROCEDURES                       !
    !--------------------------------------------------------------------------!

    pure subroutine spin_assr_rscalar(this, other)
    !!  Imports data to a spin object from a real scalar.
        class(spin), intent(inout) :: this
        real(wp),    intent(in)    :: other

        this%matrix = other*pauli0%matrix
    end subroutine

    pure subroutine spin_assr_cscalar(this, other)
    !!  Imports data to a spin object from a complex scalar.
        class(spin), intent(inout) :: this
        complex(wp), intent(in)    :: other

        this%matrix = other*pauli0%matrix
    end subroutine

    pure subroutine spin_assr_cmatrix(this, other)
    !!  Imports data to a spin object from a complex matrix.
        class(spin), intent(inout) :: this
        complex(wp), intent(in)    :: other(2, 2)

        this%matrix = other
    end subroutine

    pure subroutine spin_assr_rvector(this, other)
    !!  Imports data to a spin object from a real vector.
        class(spin), intent(inout) :: this
        real(wp),    intent(in)    :: other(8)

        ! TODO: Rewrite this without using `cx` in the future.
        this%matrix = cx(reshape(other(1:7:2), [2, 2], order=[2, 1]), &
                         reshape(other(2:8:2), [2, 2], order=[2, 1]))
    end subroutine

    !--------------------------------------------------------------------------!
    !                         SPECIFIC EXPORT PROCEDURES                       !
    !--------------------------------------------------------------------------!

    pure subroutine spin_assl_cmatrix(other, this)
    !!  Exports data from a spin object to a complex matrix.
        class(spin), intent(in)  :: this
        complex(wp), intent(out) :: other(2, 2)

        other = this%matrix
    end subroutine

    pure subroutine spin_assl_rvector(other, this)
    !!  Exports data from a spin object to a real vector.
        class(spin), intent(in)  :: this
        real(wp),    intent(out) :: other(8)

        other(1:7:2) = re([this%matrix(1, :), this%matrix(2, :)])
        other(2:8:2) = im([this%matrix(1, :), this%matrix(2, :)])
    end subroutine

    !--------------------------------------------------------------------------!
    !                     SPECIFIC EXPONENTIATION PROCEDURES                   !
    !--------------------------------------------------------------------------!

    pure function spin_expr_iscalar(this, other) result(r)
    !!  Exponentiates the spin object, where the power is a positive integer.
        class(spin), intent(in) :: this
        integer,     intent(in) :: other
        type(spin)              :: r

        integer :: n

        r = this
        do n = 2, other
            r%matrix = r%matrix*this
        end do
    end function

    !--------------------------------------------------------------------------!
    !                     SPECIFIC MULTIPLICATION PROCEDURES                   !
    !--------------------------------------------------------------------------!

    elemental function spin_mul_spin(this, other) result(r)
    !!  Defines multiplication of two spin matrices.
        class(spin), intent(in) :: this
        class(spin), intent(in) :: other
        type(spin)              :: r

        r%matrix = matmul(this%matrix, other%matrix)
    end function

    pure function spin_mull_rscalar(other, this) result(r)
    !!  Defines left multiplication of a spin matrix by a real scalar.
        class(spin), intent(in) :: this
        real(wp),    intent(in) :: other
        type(spin)              :: r

        r%matrix = other*this%matrix
    end function

    pure function spin_mulr_rscalar(this, other) result(r)
    !!  Defines right multiplication of a spin matrix by a real scalar.
        class(spin), intent(in) :: this
        real(wp),    intent(in) :: other
        type(spin)              :: r

        r%matrix = this%matrix*other
    end function

    pure function spin_mull_cscalar(other, this) result(r)
    !!  Defines left multiplication of a spin matrix by a complex scalar.
        class(spin), intent(in) :: this
        complex(wp), intent(in) :: other
        type(spin)              :: r

        r%matrix = other*this%matrix
    end function

    function spin_mulr_cscalar(this, other) result(r)
    !!  Defines right multiplication of a spin matrix by a complex scalar.
        class(spin), intent(in) :: this
        complex(wp), intent(in) :: other
        type(spin)              :: r

        r%matrix = this%matrix*other
    end function

    pure function spin_mull_cmatrix(other, this) result(r)
    !!  Defines left multiplication of a spin matrix by a complex matrix.
        class(spin), intent(in) :: this
        complex(wp), intent(in) :: other(2, 2)
        type(spin)              :: r

        r%matrix = matmul(other, this%matrix)
    end function

    pure function spin_mulr_cmatrix(this, other) result(r)
    !!  Defines right multiplication of a spin matrix by a complex matrix.
        class(spin), intent(in) :: this
        complex(wp), intent(in) :: other(2, 2)
        type(spin)              :: r

        r%matrix = matmul(this%matrix, other)
    end function

    !--------------------------------------------------------------------------!
    !                        SPECIFIC DIVISION PROCEDURES                      !
    !--------------------------------------------------------------------------!

    pure function spin_divr_rscalar(this, other) result(r)
    !!  Defines division of a spin matrix by a real scalar.
        class(spin), intent(in) :: this
        real(wp),    intent(in) :: other
        type(spin)              :: r

        r%matrix = this%matrix/other
    end function

    pure function spin_divr_cscalar(this, other) result(r)
    !!  Defines division of a spin matrix by a complex scalar.
        class(spin), intent(in) :: this
        complex(wp), intent(in) :: other
        type(spin)              :: r

        r%matrix = this%matrix/other
    end function

    !--------------------------------------------------------------------------!
    !                        SPECIFIC ADDITION PROCEDURES                      !
    !--------------------------------------------------------------------------!

    elemental function spin_add_spin(this, other) result(r)
    !!  Defines addition of two spin matrices.
        class(spin), intent(in) :: this
        class(spin), intent(in) :: other
        type(spin)              :: r

        r%matrix = this%matrix + other%matrix
    end function

    pure function spin_addl_rscalar(other, this) result(r)
    !!  Defines left addition of a spin matrix and a real scalar.
        class(spin), intent(in) :: this
        real(wp),    intent(in) :: other
        type(spin)              :: r

        r%matrix = other*pauli0%matrix + this%matrix
    end function

    pure function spin_addr_rscalar(this, other) result(r)
    !!  Defines right addition of a spin matrix and a real scalar.
        class(spin), intent(in) :: this
        real(wp),    intent(in) :: other
        type(spin)              :: r

        r%matrix = this%matrix + other*pauli0%matrix
    end function

    pure function spin_addl_cscalar(other, this) result(r)
    !!  Defines left addition of a spin matrix and a complex scalar.
        class(spin), intent(in) :: this
        complex(wp), intent(in) :: other
        type(spin)              :: r

        r%matrix = other*pauli0%matrix + this%matrix
    end function

    pure function spin_addr_cscalar(this, other) result(r)
    !!  Defines right addition of a spin matrix and a complex scalar.
        class(spin), intent(in) :: this
        complex(wp), intent(in) :: other
        type(spin)              :: r

        r%matrix = this%matrix + other*pauli0%matrix
    end function

    pure function spin_addl_cmatrix(other, this) result(r)
    !!  Defines left addition of a spin matrix and a complex matrix.
        class(spin), intent(in) :: this
        complex(wp), intent(in) :: other(2, 2)
        type(spin)              :: r

        r%matrix = other + this%matrix
    end function

    pure function spin_addr_cmatrix(this, other) result(r)
    !!  Defines right addition of a spin matrix and a complex matrix.
        class(spin), intent(in) :: this
        complex(wp), intent(in) :: other(2, 2)
        type(spin)              :: r

        r%matrix = this%matrix + other
    end function

    !--------------------------------------------------------------------------!
    !                       SPECIFIC SUBTRACTION PROCEDURES                    !
    !--------------------------------------------------------------------------!

    elemental function spin_sub_spin(this, other) result(r)
    !!  Defines subtraction of two spin matrices.
        class(spin), intent(in) :: this
        class(spin), intent(in) :: other
        type(spin)              :: r

        r%matrix = this%matrix - other%matrix
    end function

    pure function spin_subl_rscalar(other, this) result(r)
    !!  Defines left subtraction of a spin matrix and a real scalar.
        class(spin), intent(in) :: this
        real(wp),    intent(in) :: other
        type(spin)              :: r

        r%matrix = other*pauli0%matrix - this%matrix
    end function

    pure function spin_subr_rscalar(this, other) result(r)
    !!  Defines right subtraction of a spin matrix and a real scalar.
        class(spin), intent(in) :: this
        real(wp),    intent(in) :: other
        type(spin)              :: r

        r%matrix = this%matrix - other*pauli0%matrix
    end function

    pure function spin_subl_cscalar(other, this) result(r)
    !!  Defines left subtraction of a spin matrix and a complex scalar.
        class(spin), intent(in) :: this
        complex(wp), intent(in) :: other
        type(spin)              :: r

        r%matrix = other*pauli0%matrix - this%matrix
    end function

    pure function spin_subr_cscalar(this, other) result(r)
    !!  Defines right subtraction of a spin matrix and a complex scalar.
        class(spin), intent(in) :: this
        complex(wp), intent(in) :: other
        type(spin)              :: r

        r%matrix = this%matrix - other*pauli0%matrix
    end function

    pure function spin_subl_cmatrix(other, this) result(r)
    !!  Defines left subtraction of a spin matrix and a complex matrix.
        class(spin), intent(in) :: this
        complex(wp), intent(in) :: other(2, 2)
        type(spin)              :: r

        r%matrix = other - this%matrix
    end function

    pure function spin_subr_cmatrix(this, other) result(r)
    !!  Defines right subtraction of a spin matrix and a complex matrix.
        class(spin), intent(in) :: this
        complex(wp), intent(in) :: other(2, 2)
        type(spin)              :: r

        r%matrix = this%matrix - other
    end function

    !--------------------------------------------------------------------------!
    !                                MATRIX ALGEBRA                            !
    !--------------------------------------------------------------------------!

    elemental function spin_norm(this) result(r)
    !!  Calculate the Frobenius norm of the spin matrix.
        class(spin), intent(in) :: this
        real(wp)                :: r, w(8)

        w = this
        r = norm2(w)
    end function

    elemental function spin_conjg(this) result(r)
    !!  Calculate the complex conjugate of the spin matrix.
        class(spin), intent(in) :: this
        type(spin)              :: r

        r%matrix = conjg(this%matrix)
    end function

    elemental function spin_trace(this) result(r)
    !!  Calculate the trace of the spin matrix.
        class(spin), intent(in) :: this
        complex(wp)             :: r

        r = this%matrix(1, 1) + this%matrix(2, 2)
    end function

    pure function spin_inv(this) result(r)
    !!  Calculate the inverse of the spin matrix.
        class(spin), intent(in) :: this
        type(spin)              :: r

        r%matrix = inverse(this%matrix)
    end function

    pure function spin_sum(this) result(r)
    !!  Calculate the sum of an array of spin matrices.
        class(spin), intent(in) :: this(:)
        type(spin)              :: r
        integer                 :: n

        do n = 1, size(this)
            r%matrix = r%matrix + this(n)%matrix
        end do
    end function
end module
