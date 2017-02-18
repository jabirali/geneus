!> Author:   Jabir Ali Ouassou
!> Date:     2015-07-10
!> Category: Foundation
!>
!> This module defines the data type 'nambu', which represents 4×4 complex matrices in spin and Nambu space. The module
!> overloads arithmetic operators to work with the new type, and defines and exports the Pauli matrices as constants.

module nambu_m
  use :: math_m
  use :: spin_m
  use :: matrix_m
  private

  ! Public interface
  public nambu, nambuv
  public inverse, trace, conjg, sum

  ! Type declaration
  type nambu
    complex(wp) :: matrix(4,4) = 0.0_wp   !! Encapsulate a spin and Nambu space matrix
  contains
    ! Overload constructors and operators
    generic :: nambu         => cons_rscalar, cons_cscalar, cons_cmatrix, cons_nambu
    generic :: assignment(=) => assr_rscalar, assr_cscalar, assr_cmatrix, assl_cmatrix
    generic :: operator(+)   => addl_rscalar, addr_rscalar, addl_cscalar, addr_cscalar, addl_cmatrix, addr_cmatrix, add_nambu
    generic :: operator(-)   => subl_rscalar, subr_rscalar, subl_cscalar, subr_cscalar, subl_cmatrix, subr_cmatrix, sub_nambu
    generic :: operator(*)   => mull_rscalar, mulr_rscalar, mull_cscalar, mulr_cscalar, mull_cmatrix, mulr_cmatrix, mul_nambu
    generic :: operator(/)   => divr_rscalar, divr_cscalar
    generic :: operator(**)  => expr_iscalar

    ! Specific methods for construction
    procedure,   nopass,   private :: cons_nambu   => nambu_cons_nambu     !! Construction from a nambu object
    procedure,   nopass,   private :: cons_rscalar => nambu_cons_rscalar   !! Construction from a real scalar 
    procedure,   nopass,   private :: cons_cscalar => nambu_cons_cscalar   !! Construction from a complex scalar 
    procedure,   nopass,   private :: cons_cmatrix => nambu_cons_cmatrix   !! Construction from a complex matrix

    ! Specific implementations of assignments
    procedure, pass(this), private :: assr_rscalar => nambu_assr_rscalar   !! Assign data from a real scalar      (right)
    procedure, pass(this), private :: assr_cscalar => nambu_assr_cscalar   !! Assign data from a complex scalar   (right)
    procedure, pass(this), private :: assl_cmatrix => nambu_assl_cmatrix   !! Assign data into a complex matrix   (left)
    procedure, pass(this), private :: assr_cmatrix => nambu_assr_cmatrix   !! Assign data from a complex matrix   (right)

    ! Specific implementations of addition
    procedure, pass(this), private :: add_nambu    => nambu_add_nambu      !! Addition with a nambu object
    procedure, pass(this), private :: addl_rscalar => nambu_addl_rscalar   !! Addition with a real scalar         (left) 
    procedure, pass(this), private :: addr_rscalar => nambu_addr_rscalar   !! Addition with a real scalar         (right)
    procedure, pass(this), private :: addl_cscalar => nambu_addl_cscalar   !! Addition with a complex scalar      (left) 
    procedure, pass(this), private :: addr_cscalar => nambu_addr_cscalar   !! Addition with a complex scalar      (right)
    procedure, pass(this), private :: addl_cmatrix => nambu_addl_cmatrix   !! Addition with a complex matrix      (left) 
    procedure, pass(this), private :: addr_cmatrix => nambu_addr_cmatrix   !! Addition with a complex matrix      (right)

    ! Specific implementations of subtraction
    procedure, pass(this), private :: sub_nambu    => nambu_sub_nambu      !! Subtraction with a nambu object
    procedure, pass(this), private :: subl_rscalar => nambu_subl_rscalar   !! Subtraction with a real scalar      (left) 
    procedure, pass(this), private :: subr_rscalar => nambu_subr_rscalar   !! Subtraction with a real scalar      (right)
    procedure, pass(this), private :: subl_cscalar => nambu_subl_cscalar   !! Subtraction with a complex scalar   (left)
    procedure, pass(this), private :: subr_cscalar => nambu_subr_cscalar   !! Subtraction with a complex scalar   (right) 
    procedure, pass(this), private :: subl_cmatrix => nambu_subl_cmatrix   !! Subtraction with a complex matrix   (left) 
    procedure, pass(this), private :: subr_cmatrix => nambu_subr_cmatrix   !! Subtraction with a complex matrix   (right) 

    ! Specific implementations of multiplication
    procedure, pass(this), private :: mul_nambu    => nambu_mul_nambu      !! Multiplication by a nambu object
    procedure, pass(this), private :: mull_rscalar => nambu_mull_rscalar   !! Multiplication by a real scalar     (left) 
    procedure, pass(this), private :: mulr_rscalar => nambu_mulr_rscalar   !! Multiplication by a real scalar     (right)
    procedure, pass(this), private :: mull_cscalar => nambu_mull_cscalar   !! Multiplication by a complex scalar  (left)
    procedure, pass(this), private :: mulr_cscalar => nambu_mulr_cscalar   !! Multiplication by a complex scalar  (right)
    procedure, pass(this), private :: mull_cmatrix => nambu_mull_cmatrix   !! Multiplication by a complex matrix  (left)
    procedure, pass(this), private :: mulr_cmatrix => nambu_mulr_cmatrix   !! Multiplication by a complex matrix  (right)

    ! Specific implementations of division
    procedure, pass(this), private :: divr_rscalar => nambu_divr_rscalar   !! Division by a real scalar           (right)
    procedure, pass(this), private :: divr_cscalar => nambu_divr_cscalar   !! Division by a complex scalar        (right)

    ! Specific implementations of exponentiation
    procedure, pass(this), private :: expr_iscalar => nambu_expr_iscalar   !! Exponentiation by an integer        (right)
  end type

  ! Matrix inverse
  interface inverse
    module procedure nambu_inv
  end interface

  ! Matrix trace
  interface trace
    module procedure nambu_trace
  end interface

  ! Matrix sums
  interface sum
    module procedure nambu_sum
  end interface

  ! Complex conjugation
  interface conjg
    module procedure nambu_conjg
  end interface
contains

  !--------------------------------------------------------------------------------!
  !                            SPECIFIC CONSTRUCTORS                               !
  !--------------------------------------------------------------------------------!

  pure function nambuv(n) result(r)
    !! Constructs basis vector number n in spin-nambu space.
    integer, intent(in) :: n
    type(nambu)         :: r

    select case (n)
      case(0)
        ! Basis matrix τ₀σ₀
        r % matrix(1:2,1:2) = +pauli0 % matrix
        r % matrix(3:4,3:4) = +pauli0 % matrix
      case(1)
        ! Basis matrix τ₀σ₁
        r % matrix(1:2,1:2) = +pauli1 % matrix
        r % matrix(3:4,3:4) = +pauli1 % matrix
      case(2)
        ! Basis matrix τ₀σ₂
        r % matrix(1:2,1:2) = +pauli2 % matrix
        r % matrix(3:4,3:4) = -pauli2 % matrix
      case(3)
        ! Basis matrix τ₀σ₃
        r % matrix(1:2,1:2) = +pauli3 % matrix
        r % matrix(3:4,3:4) = +pauli3 % matrix
      case(4)
        ! Basis matrix τ₃σ₀
        r % matrix(1:2,1:2) = +pauli0 % matrix
        r % matrix(3:4,3:4) = -pauli0 % matrix
      case(5)
        ! Basis matrix τ₃σ₁
        r % matrix(1:2,1:2) = +pauli1 % matrix
        r % matrix(3:4,3:4) = -pauli1 % matrix
      case(6)
        ! Basis matrix τ₃σ₂
        r % matrix(1:2,1:2) = +pauli2 % matrix
        r % matrix(3:4,3:4) = +pauli2 % matrix
      case(7)
        ! Basis matrix τ₃σ₃
        r % matrix(1:2,1:2) = +pauli3 % matrix
        r % matrix(3:4,3:4) = -pauli3 % matrix
    end select
  end function

  pure function nambu_cons_rscalar(other) result(this)
    !! Constructs a nambu object from a real scalar.
    real(wp),   intent(in) :: other
    type(nambu)            :: this

    this = other
  end function

  pure function nambu_cons_cscalar(other) result(this)
    !! Constructs a nambu object from a complex scalar.
    complex(wp), intent(in)  :: other
    type(nambu)              :: this

    this = other
  end function

  pure function nambu_cons_cmatrix(other) result(this)
    !! Constructs a nambu object from a complex matrix.
    complex(wp), intent(in) :: other(4,4)
    type(nambu)             :: this

    this = other
  end function

  pure function nambu_cons_nambu(other) result(this)
    !! Constructs a nambu object from an existing one.
    type(nambu), intent(in) :: other
    type(nambu)             :: this

    this = other
  end function

  !--------------------------------------------------------------------------------!
  !                         SPECIFIC IMPORT PROCEDURES                             !
  !--------------------------------------------------------------------------------!

  pure subroutine nambu_assr_rscalar(this, other)
    !! Imports data to a nambu object from a real scalar.
    class(nambu), intent(inout) :: this
    real(wp),     intent(in)    :: other

    this%matrix = other * identity4
  end subroutine

  pure subroutine nambu_assr_cscalar(this, other)
    !! Imports data to a nambu object from a complex scalar.
    class(nambu), intent(inout) :: this
    complex(wp),  intent(in)    :: other

    this%matrix = other * identity4
  end subroutine

  pure subroutine nambu_assr_cmatrix(this, other)
    !! Imports data to a nambu object from a complex matrix.
    class(nambu),  intent(inout) :: this
    complex(wp),   intent(in)    :: other(4,4)

    this % matrix = other 
  end subroutine



  !--------------------------------------------------------------------------------!
  !                         SPECIFIC EXPORT PROCEDURES                             !
  !--------------------------------------------------------------------------------!

  pure subroutine nambu_assl_cmatrix(other, this)
    !! Exports data from a nambu object to a complex matrix.
    class(nambu), intent(in)  :: this
    complex(wp),  intent(out) :: other(4,4)

    other = this%matrix 
  end subroutine



  !--------------------------------------------------------------------------------!
  !                     SPECIFIC EXPONENTIATION PROCEDURES                         !
  !--------------------------------------------------------------------------------!

  pure function nambu_expr_iscalar(this, other) result(r)
    !! Exponentiates the nambu object, where the power is a positive integer.
    class(nambu), intent(in) :: this
    integer,      intent(in) :: other
    type(nambu)              :: r
    integer                  :: n

    r = this
    do n=2,other
      r%matrix = r%matrix * this
    end do
  end function



  !--------------------------------------------------------------------------------!
  !                     SPECIFIC MULTIPLICATION PROCEDURES                         !
  !--------------------------------------------------------------------------------!

  elemental pure function nambu_mul_nambu(this, other) result(r)
    !! Defines multiplication of two nambu matrices.
    class(nambu), intent(in) :: this
    class(nambu), intent(in) :: other
    type(nambu)              :: r

    r%matrix = matmul(this%matrix, other%matrix)
  end function

  pure function nambu_mull_rscalar(other, this) result(r)
    !! Defines left multiplication of a nambu matrix by a real scalar.
    class(nambu), intent(in) :: this
    real(wp),     intent(in) :: other
    type(nambu)              :: r

    r%matrix = other * this%matrix
  end function

  pure function nambu_mulr_rscalar(this, other) result(r)
    !! Defines right multiplication of a nambu matrix by a real scalar.
    class(nambu), intent(in) :: this
    real(wp),     intent(in) :: other
    type(nambu)              :: r

    r%matrix = this%matrix * other
  end function

  pure function nambu_mull_cscalar(other, this) result(r)
    !! Defines left multiplication of a nambu matrix by a complex scalar.
    class(nambu), intent(in) :: this
    complex(wp),  intent(in) :: other
    type(nambu)              :: r

    r%matrix = other * this%matrix
  end function

  function nambu_mulr_cscalar(this, other) result(r)
    !! Defines right multiplication of a nambu matrix by a complex scalar.
    class(nambu), intent(in) :: this
    complex(wp),  intent(in) :: other
    type(nambu)              :: r

    r%matrix = this%matrix * other
  end function

  pure function nambu_mull_cmatrix(other, this) result(r)
    !! Defines left multiplication of a nambu matrix by a complex matrix.
    class(nambu), intent(in) :: this
    complex(wp),  intent(in) :: other(4,4)
    type(nambu)              :: r
 
    r%matrix = matmul(other, this%matrix)
  end function

  pure function nambu_mulr_cmatrix(this, other) result(r)
    !! Defines right multiplication of a nambu matrix by a complex matrix.
    class(nambu), intent(in) :: this
    complex(wp),  intent(in) :: other(4,4)
    type(nambu)              :: r

    r%matrix = matmul(this%matrix, other)
  end function



  !--------------------------------------------------------------------------------!
  !                        SPECIFIC DIVISION PROCEDURES                            !
  !--------------------------------------------------------------------------------!

  pure function nambu_divr_rscalar(this, other) result(r)
    !! Defines division of a nambu matrix by a real scalar.
    class(nambu), intent(in) :: this
    real(wp),     intent(in) :: other
    type(nambu)              :: r

    r%matrix = this%matrix / other
  end function

  pure function nambu_divr_cscalar(this, other) result(r)
    !! Defines division of a nambu matrix by a complex scalar.
    class(nambu), intent(in) :: this
    complex(wp),  intent(in) :: other
    type(nambu)              :: r

    r%matrix = this%matrix / other
  end function



  !--------------------------------------------------------------------------------!
  !                        SPECIFIC ADDITION PROCEDURES                            !
  !--------------------------------------------------------------------------------!

  elemental pure function nambu_add_nambu(this, other) result(r)
    !! Defines addition of two nambu matrices.
    class(nambu), intent(in) :: this
    class(nambu), intent(in) :: other
    type(nambu)              :: r

    r%matrix = this%matrix + other%matrix
  end function

  pure function nambu_addl_rscalar(other, this) result(r)
    !! Defines left addition of a nambu matrix and a real scalar.
    class(nambu), intent(in) :: this
    real(wp),     intent(in) :: other
    type(nambu)              :: r

    r%matrix = other*identity4 + this%matrix
  end function

  pure function nambu_addr_rscalar(this, other) result(r)
    !! Defines right addition of a nambu matrix and a real scalar.
    class(nambu), intent(in) :: this
    real(wp),     intent(in) :: other
    type(nambu)              :: r

    r%matrix = this%matrix + other*identity4
  end function

  pure function nambu_addl_cscalar(other, this) result(r)
    !! Defines left addition of a nambu matrix and a complex scalar.
    class(nambu), intent(in) :: this
    complex(wp),  intent(in) :: other
    type(nambu)              :: r

    r%matrix = other*identity4 + this%matrix
  end function

  pure function nambu_addr_cscalar(this, other) result(r)
    !! Defines right addition of a nambu matrix and a complex scalar.
    class(nambu), intent(in) :: this
    complex(wp),  intent(in) :: other
    type(nambu)              :: r

    r%matrix = this%matrix + other*identity4
  end function

  pure function nambu_addl_cmatrix(other, this) result(r)
    !! Defines left addition of a nambu matrix and a complex matrix.
    class(nambu), intent(in) :: this
    complex(wp),  intent(in) :: other(4,4)
    type(nambu)              :: r

    r%matrix = other + this%matrix
  end function

  pure function nambu_addr_cmatrix(this, other) result(r)
    !! Defines right addition of a nambu matrix and a complex matrix.
    class(nambu), intent(in) :: this
    complex(wp),  intent(in) :: other(4,4)
    type(nambu)              :: r

    r%matrix = this%matrix + other
  end function



  !--------------------------------------------------------------------------------!
  !                       SPECIFIC SUBTRACTION PROCEDURES                          !
  !--------------------------------------------------------------------------------!

  elemental pure function nambu_sub_nambu(this, other) result(r)
    !! Defines subtraction of two nambu matrices.
    class(nambu), intent(in) :: this
    class(nambu), intent(in) :: other
    type(nambu)              :: r

    r%matrix = this%matrix - other%matrix
  end function

  pure function nambu_subl_rscalar(other, this) result(r)
    !! Defines left subtraction of a nambu matrix and a real scalar.
    class(nambu), intent(in) :: this
    real(wp),     intent(in) :: other
    type(nambu)              :: r

    r%matrix = other*identity4 - this%matrix
  end function

  pure function nambu_subr_rscalar(this, other) result(r)
    !! Defines right subtraction of a nambu matrix and a real scalar.
    class(nambu), intent(in) :: this
    real(wp),     intent(in) :: other
    type(nambu)              :: r

    r%matrix = this%matrix - other*identity4
  end function

  pure function nambu_subl_cscalar(other, this) result(r)
    !! Defines left subtraction of a nambu matrix and a complex scalar.
    class(nambu), intent(in) :: this
    complex(wp),  intent(in) :: other
    type(nambu)              :: r

    r%matrix = other*identity4 - this%matrix
  end function

  pure function nambu_subr_cscalar(this, other) result(r)
    !! Defines right subtraction of a nambu matrix and a complex scalar.
    class(nambu), intent(in) :: this
    complex(wp),  intent(in) :: other
    type(nambu)              :: r

    r%matrix = this%matrix - other*identity4
  end function

  pure function nambu_subl_cmatrix(other, this) result(r)
    !! Defines left subtraction of a nambu matrix and a complex matrix.
    class(nambu), intent(in) :: this
    complex(wp),  intent(in) :: other(4,4)
    type(nambu)              :: r

    r%matrix = other - this%matrix
  end function

  pure function nambu_subr_cmatrix(this, other) result(r)
    !! Defines right subtraction of a nambu matrix and a complex matrix.
    class(nambu), intent(in) :: this
    complex(wp),  intent(in) :: other(4,4)
    type(nambu)              :: r

    r%matrix = this%matrix - other
  end function



  !--------------------------------------------------------------------------------!
  !                                MATRIX ALGEBRA                                  !
  !--------------------------------------------------------------------------------!

  elemental pure function nambu_conjg(this) result(r)
    !! Calculate the complex conjugate of the nambu matrix.
    class(nambu), intent(in) :: this
    type(nambu)              :: r

    r%matrix = conjg(this%matrix)
  end function

  elemental pure function nambu_trace(this) result(r)
    !! Calculate the trace of the nambu matrix.
    class(nambu), intent(in) :: this
    complex(wp)              :: r

    r = this%matrix(1,1) + this%matrix(2,2) + this%matrix(3,3) + this%matrix(4,4)
  end function

  pure function nambu_inv(this) result(r)
    !! Calculate the inverse of the nambu matrix.
    class(nambu), intent(in) :: this
    type(nambu)              :: r

    r%matrix = inverse(this%matrix)
  end function

  pure function nambu_sum(this) result(r)
    !! Calculate the sum of an array of nambu matrices.
    class(nambu), intent(in) :: this(:)
    type(nambu)              :: r
    integer                 :: n

    do n=1,size(this)
      r % matrix = r % matrix + this(n) % matrix
    end do
  end function
end module
