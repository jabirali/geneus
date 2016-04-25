!> Author:   Jabir Ali Ouassou
!> Date:     2015-07-10
!> Category: Foundation
!>
!> This module defines the data type 'spin', which represents 2×2 complex matrices in spin space. The module overloads
!> common arithmetic operators to work with the new data type, and defines and exports the Pauli matrices as constants.
!> To make it easier to interact with common differential equation solvers, which often operate on real state vectors,
!> the assignment operator is overloaded in such a way that 'spin' can be easily imported/exported to a real vector(8).

module spin_m
  use :: math_m
  private

  ! Public interface
  public spin
  public inv, trace, conjg, norm2, sum
  public pauli, pauli0, pauli1, pauli2, pauli3

  ! Type declaration
  type spin
    complex(wp)     :: matrix(2,2)   =  0.0_wp
  contains
    ! Overload constructor
    generic, public :: spin          => construct_rscalar_, &
                                        construct_cscalar_, &
                                        construct_cmatrix_, &
                                        construct_rvector_, &
                                        construct_spin_ 
    ! Overload assignment
    generic, public :: assignment(=) => import_rscalar_,    &
                                        import_cscalar_,    &
                                        import_cmatrix_,    &
                                        import_rvector_,    &
                                        export_cmatrix_,    &
                                        export_rvector_ 

    ! Overload exponentiation
    generic, public :: operator(**)  => exp_integer_

    ! Overload multiplication
    generic, public :: operator(*)   => multl_rscalar_,     &
                                        multr_rscalar_,     &
                                        multl_cscalar_,     &
                                        multr_cscalar_,     &
                                        multl_cmatrix_,     &
                                        multr_cmatrix_,     &
                                        mult_spin_ 

    ! Overload division
    generic, public :: operator(/)   => div_rscalar_,       &
                                        div_cscalar_
    ! Overload addition
    generic, public :: operator(+)   => addl_rscalar_,      &
                                        addr_rscalar_,      &
                                        addl_cscalar_,      &
                                        addr_cscalar_,      &
                                        addl_cmatrix_,      &
                                        addr_cmatrix_,      &
                                        add_spin_ 
    ! Overload subtraction
    generic, public :: operator(-)   => subl_rscalar_,      &
                                        subr_rscalar_,      &
                                        subl_cscalar_,      &
                                        subr_cscalar_,      &
                                        subl_cmatrix_,      &
                                        subr_cmatrix_,      &
                                        sub_spin_ 

    ! Specific methods for construction
    procedure, nopass,     private :: construct_spin_    => spin_construct_spin      !! Construction from a spin object
    procedure, nopass,     private :: construct_rscalar_ => spin_construct_rscalar   !! Construction from a real scalar 
    procedure, nopass,     private :: construct_cscalar_ => spin_construct_cscalar   !! Construction from a complex scalar 
    procedure, nopass,     private :: construct_cmatrix_ => spin_construct_cmatrix   !! Construction from a complex matrix
    procedure, nopass,     private :: construct_rvector_ => spin_construct_rvector   !! Construction from a real vector

    ! Specific methods for importing the object state
    procedure, pass(this), private :: import_rscalar_    => spin_import_rscalar      !! Import data from a real scalar 
    procedure, pass(this), private :: import_cscalar_    => spin_import_cscalar      !! Import data from a complex scalar
    procedure, pass(this), private :: import_cmatrix_    => spin_import_cmatrix      !! Import data from a complex matrix
    procedure, pass(this), private :: import_rvector_    => spin_import_rvector      !! Import data from a real vector 
 
    ! Private methods for exporting the object state 
    procedure, pass(this), private :: export_cmatrix_    => spin_export_cmatrix      !! Export data to a complex matrix
    procedure, pass(this), private :: export_rvector_    => spin_export_rvector      !! Export data to a real vector

    ! Specific implementations of exponentiation
    procedure, pass(this), private :: exp_integer_       => spin_exp_integer         !! Exponentiation by an integer

    ! Specific implementations of multiplication
    procedure, pass(this), private :: mult_spin_         => spin_mult_spin           !! Multiplication by a spin object
    procedure, pass(this), private :: multl_rscalar_     => spin_multl_rscalar       !! Multiplication by a real scalar (left) 
    procedure, pass(this), private :: multr_rscalar_     => spin_multr_rscalar       !! Multiplication by a real scalar (right)
    procedure, pass(this), private :: multl_cscalar_     => spin_multl_cscalar       !! Multiplication by a complex scalar (left)
    procedure, pass(this), private :: multr_cscalar_     => spin_multr_cscalar       !! Multiplication by a complex scalar (right)
    procedure, pass(this), private :: multl_cmatrix_     => spin_multl_cmatrix       !! Multiplication by a complex matrix (left)
    procedure, pass(this), private :: multr_cmatrix_     => spin_multr_cmatrix       !! Multiplication by a complex matrix (right)

    ! Specific implementations of division
    procedure, pass(this), private :: div_rscalar_       => spin_div_rscalar        !! Division by a real scalar
    procedure, pass(this), private :: div_cscalar_       => spin_div_cscalar        !! Division by a complex scalar

    ! Specific implementations of addition
    procedure, pass(this), private :: add_spin_          => spin_add_spin           !! Addition with a spin object
    procedure, pass(this), private :: addl_rscalar_      => spin_addl_rscalar       !! Addition with a real scalar (left) 
    procedure, pass(this), private :: addr_rscalar_      => spin_addr_rscalar       !! Addition with a real scalar (right)
    procedure, pass(this), private :: addl_cscalar_      => spin_addl_cscalar       !! Addition with a complex scalar (left) 
    procedure, pass(this), private :: addr_cscalar_      => spin_addr_cscalar       !! Addition with a complex scalar (right)
    procedure, pass(this), private :: addl_cmatrix_      => spin_addl_cmatrix       !! Addition with a complex matrix (left) 
    procedure, pass(this), private :: addr_cmatrix_      => spin_addr_cmatrix       !! Addition with a complex matrix (right)

    ! Specific implementations of subtraction
    procedure, pass(this), private :: sub_spin_          => spin_sub_spin           !! Subtraction with a spin object
    procedure, pass(this), private :: subl_rscalar_      => spin_subl_rscalar       !! Subtraction with a real scalar (left) 
    procedure, pass(this), private :: subr_rscalar_      => spin_subr_rscalar       !! Subtraction with a real scalar (right)
    procedure, pass(this), private :: subl_cscalar_      => spin_subl_cscalar       !! Subtraction with a complex scalar (left)
    procedure, pass(this), private :: subr_cscalar_      => spin_subr_cscalar       !! Subtraction with a complex scalar (right) 
    procedure, pass(this), private :: subl_cmatrix_      => spin_subl_cmatrix       !! Subtraction with a complex matrix (left) 
    procedure, pass(this), private :: subr_cmatrix_      => spin_subr_cmatrix       !! Subtraction with a complex matrix (right) 
  end type

  ! Matrix inverse
  interface inv
    module procedure spin_inv
  end interface

  ! Matrix trace
  interface trace
    module procedure spin_trace
  end interface

  ! Matrix sums
  interface sum
    module procedure spin_sum
  end interface

  ! Complex conjugation
  interface conjg
    module procedure spin_conjg
  end interface

  ! Matrix norm
  interface norm2
    module procedure spin_norm
  end interface

  ! Exported constants
  type(spin), parameter :: pauli0     = spin(reshape([ ( 1, 0), ( 0, 0), ( 0, 0), ( 1, 0) ], [2,2], order=[2,1]))
  type(spin), parameter :: pauli1     = spin(reshape([ ( 0, 0), ( 1, 0), ( 1, 0), ( 0, 0) ], [2,2], order=[2,1]))
  type(spin), parameter :: pauli2     = spin(reshape([ ( 0, 0), ( 0,-1), ( 0, 1), ( 0, 0) ], [2,2], order=[2,1]))
  type(spin), parameter :: pauli3     = spin(reshape([ ( 1, 0), ( 0, 0), ( 0, 0), (-1, 0) ], [2,2], order=[2,1]))
  type(spin), parameter :: pauli(0:3) = [pauli0, pauli1, pauli2, pauli3]
contains

  !--------------------------------------------------------------------------------!
  !                            SPECIFIC CONSTRUCTORS                               !
  !--------------------------------------------------------------------------------!

  pure function spin_construct_rscalar(other) result(this)
    !! Constructs a spin object from a real scalar.
    real(wp),  intent(in)  :: other
    type(spin)             :: this

    this = other
  end function

  pure function spin_construct_cscalar(other) result(this)
    !! Constructs a spin object from a complex scalar.
    complex(wp), intent(in)  :: other
    type(spin)               :: this

    this = other
  end function

  pure function spin_construct_cmatrix(other) result(this)
    !! Constructs a spin object from a complex matrix.
    complex(wp), intent(in) :: other(2,2)
    type(spin)              :: this

    this = other
  end function

  pure function spin_construct_rvector(other) result(this)
    !! Constructs a spin object from a real vector.
    real(wp), intent(in) :: other(8)
    type(spin)           :: this

    this = other
  end function

  pure function spin_construct_spin(other) result(this)
    !! Constructs a spin object from an existing one.
    type(spin), intent(in) :: other
    type(spin)             :: this

    this = other
  end function



  !--------------------------------------------------------------------------------!
  !                         SPECIFIC IMPORT PROCEDURES                             !
  !--------------------------------------------------------------------------------!

  pure subroutine spin_import_rscalar(this, other)
    !! Imports data to a spin object from a real scalar.
    class(spin), intent(inout) :: this
    real(wp),    intent(in)    :: other

    this%matrix = other * pauli0%matrix
  end subroutine

  pure subroutine spin_import_cscalar(this, other)
    !! Imports data to a spin object from a complex scalar.
    class(spin), intent(inout) :: this
    complex(wp), intent(in)    :: other

    this%matrix = other * pauli0%matrix
  end subroutine

  pure subroutine spin_import_cmatrix(this, other)
    !! Imports data to a spin object from a complex matrix.
    class(spin),  intent(inout) :: this
    complex(wp),  intent(in)    :: other(2,2)

    this%matrix = other
  end subroutine

  pure subroutine spin_import_rvector(this, other)
    !! Imports data to a spin object from a real vector.
    class(spin), intent(inout) :: this
    real(wp),    intent(in)    :: other(8)

    this%matrix = cx(reshape(other(1:7:2),[2,2],order=[2,1]),&
                     reshape(other(2:8:2),[2,2],order=[2,1]))
  end subroutine



  !--------------------------------------------------------------------------------!
  !                         SPECIFIC EXPORT PROCEDURES                             !
  !--------------------------------------------------------------------------------!

  pure subroutine spin_export_cmatrix(other, this)
    !! Exports data from a spin object to a complex matrix.
    class(spin), intent(in)  :: this
    complex(wp), intent(out) :: other(2,2)

    other = this%matrix 
  end subroutine

  pure subroutine spin_export_rvector(other, this)
    !! Exports data from a spin object to a real vector.
    class(spin), intent(in)  :: this
    real(wp),    intent(out) :: other(8)

    other(1:7:2) = re([ this%matrix(1,:), this%matrix(2,:) ])
    other(2:8:2) = im([ this%matrix(1,:), this%matrix(2,:) ])
  end subroutine



  !--------------------------------------------------------------------------------!
  !                     SPECIFIC EXPONENTIATION PROCEDURES                         !
  !--------------------------------------------------------------------------------!

  pure function spin_exp_integer(this, other) result(r)
    !! Exponentiates the spin object, where the power is a positive integer.
    class(spin), intent(in) :: this
    integer,     intent(in) :: other
    type(spin)              :: r
    integer                 :: n

    r = this
    do n=2,other
      r%matrix = r%matrix * this
    end do
  end function



  !--------------------------------------------------------------------------------!
  !                     SPECIFIC MULTIPLICATION PROCEDURES                         !
  !--------------------------------------------------------------------------------!

  elemental pure function spin_mult_spin(this, other) result(r)
    !! Defines multiplication of two spin matrices.
    class(spin), intent(in) :: this
    class(spin), intent(in) :: other
    type(spin)              :: r

    r%matrix = matmul(this%matrix, other%matrix)
  end function

  pure function spin_multl_rscalar(other, this) result(r)
    !! Defines left multiplication of a spin matrix by a real scalar.
    class(spin), intent(in) :: this
    real(wp),    intent(in) :: other
    type(spin)              :: r

    r%matrix = other * this%matrix
  end function

  pure function spin_multr_rscalar(this, other) result(r)
    !! Defines right multiplication of a spin matrix by a real scalar.
    class(spin), intent(in) :: this
    real(wp),    intent(in) :: other
    type(spin)              :: r

    r%matrix = this%matrix * other
  end function

  pure function spin_multl_cscalar(other, this) result(r)
    !! Defines left multiplication of a spin matrix by a complex scalar.
    class(spin), intent(in) :: this
    complex(wp), intent(in) :: other
    type(spin)              :: r

    r%matrix = other * this%matrix
  end function

  function spin_multr_cscalar(this, other) result(r)
    !! Defines right multiplication of a spin matrix by a complex scalar.
    class(spin), intent(in) :: this
    complex(wp), intent(in) :: other
    type(spin)              :: r

    r%matrix = this%matrix * other
  end function

  pure function spin_multl_cmatrix(other, this) result(r)
    !! Defines left multiplication of a spin matrix by a complex matrix.
    class(spin), intent(in) :: this
    complex(wp), intent(in) :: other(2,2)
    type(spin)              :: r
 
    r%matrix = matmul(other, this%matrix)
  end function

  pure function spin_multr_cmatrix(this, other) result(r)
    !! Defines right multiplication of a spin matrix by a complex matrix.
    class(spin), intent(in) :: this
    complex(wp), intent(in) :: other(2,2)
    type(spin)              :: r

    r%matrix = matmul(this%matrix, other)
  end function



  !--------------------------------------------------------------------------------!
  !                        SPECIFIC DIVISION PROCEDURES                            !
  !--------------------------------------------------------------------------------!

  pure function spin_div_rscalar(this, other) result(r)
    !! Defines division of a spin matrix by a real scalar.
    class(spin), intent(in) :: this
    real(wp),    intent(in) :: other
    type(spin)              :: r

    r%matrix = this%matrix / other
  end function

  pure function spin_div_cscalar(this, other) result(r)
    !! Defines division of a spin matrix by a complex scalar.
    class(spin), intent(in) :: this
    complex(wp), intent(in) :: other
    type(spin)              :: r

    r%matrix = this%matrix / other
  end function



  !--------------------------------------------------------------------------------!
  !                        SPECIFIC ADDITION PROCEDURES                            !
  !--------------------------------------------------------------------------------!

  elemental pure function spin_add_spin(this, other) result(r)
    !! Defines addition of two spin matrices.
    class(spin), intent(in) :: this
    class(spin), intent(in) :: other
    type(spin)              :: r

    r%matrix = this%matrix + other%matrix
  end function

  pure function spin_addl_rscalar(other, this) result(r)
    !! Defines left addition of a spin matrix and a real scalar.
    class(spin), intent(in) :: this
    real(wp),    intent(in) :: other
    type(spin)              :: r

    r%matrix = other*pauli0%matrix + this%matrix
  end function

  pure function spin_addr_rscalar(this, other) result(r)
    !! Defines right addition of a spin matrix and a real scalar.
    class(spin), intent(in) :: this
    real(wp),    intent(in) :: other
    type(spin)              :: r

    r%matrix = this%matrix + other*pauli0%matrix
  end function

  pure function spin_addl_cscalar(other, this) result(r)
    !! Defines left addition of a spin matrix and a complex scalar.
    class(spin), intent(in) :: this
    complex(wp), intent(in) :: other
    type(spin)              :: r

    r%matrix = other*pauli0%matrix + this%matrix
  end function

  pure function spin_addr_cscalar(this, other) result(r)
    !! Defines right addition of a spin matrix and a complex scalar.
    class(spin), intent(in) :: this
    complex(wp), intent(in) :: other
    type(spin)              :: r

    r%matrix = this%matrix + other*pauli0%matrix
  end function

  pure function spin_addl_cmatrix(other, this) result(r)
    !! Defines left addition of a spin matrix and a complex matrix.
    class(spin), intent(in) :: this
    complex(wp), intent(in) :: other(2,2)
    type(spin)              :: r

    r%matrix = other + this%matrix
  end function

  pure function spin_addr_cmatrix(this, other) result(r)
    !! Defines right addition of a spin matrix and a complex matrix.
    class(spin), intent(in) :: this
    complex(wp), intent(in) :: other(2,2)
    type(spin)              :: r

    r%matrix = this%matrix + other
  end function



  !--------------------------------------------------------------------------------!
  !                       SPECIFIC SUBTRACTION PROCEDURES                          !
  !--------------------------------------------------------------------------------!

  elemental pure function spin_sub_spin(this, other) result(r)
    !! Defines subtraction of two spin matrices.
    class(spin), intent(in) :: this
    class(spin), intent(in) :: other
    type(spin)              :: r

    r%matrix = this%matrix - other%matrix
  end function

  pure function spin_subl_rscalar(other, this) result(r)
    !! Defines left subtraction of a spin matrix and a real scalar.
    class(spin), intent(in) :: this
    real(wp),    intent(in) :: other
    type(spin)              :: r

    r%matrix = other*pauli0%matrix - this%matrix
  end function

  pure function spin_subr_rscalar(this, other) result(r)
    !! Defines right subtraction of a spin matrix and a real scalar.
    class(spin), intent(in) :: this
    real(wp),    intent(in) :: other
    type(spin)              :: r

    r%matrix = this%matrix - other*pauli0%matrix
  end function

  pure function spin_subl_cscalar(other, this) result(r)
    !! Defines left subtraction of a spin matrix and a complex scalar.
    class(spin), intent(in) :: this
    complex(wp), intent(in) :: other
    type(spin)              :: r

    r%matrix = other*pauli0%matrix - this%matrix
  end function

  pure function spin_subr_cscalar(this, other) result(r)
    !! Defines right subtraction of a spin matrix and a complex scalar.
    class(spin), intent(in) :: this
    complex(wp), intent(in) :: other
    type(spin)              :: r

    r%matrix = this%matrix - other*pauli0%matrix
  end function

  pure function spin_subl_cmatrix(other, this) result(r)
    !! Defines left subtraction of a spin matrix and a complex matrix.
    class(spin), intent(in) :: this
    complex(wp), intent(in) :: other(2,2)
    type(spin)              :: r

    r%matrix = other - this%matrix
  end function

  pure function spin_subr_cmatrix(this, other) result(r)
    !! Defines right subtraction of a spin matrix and a complex matrix.
    class(spin), intent(in) :: this
    complex(wp), intent(in) :: other(2,2)
    type(spin)              :: r

    r%matrix = this%matrix - other
  end function



  !--------------------------------------------------------------------------------!
  !                                MATRIX ALGEBRA                                  !
  !--------------------------------------------------------------------------------!

  elemental pure function spin_norm(this) result(r)
    !! Calculate the Frobenius norm of the spin matrix.
    class(spin), intent(in) :: this
    real(wp)                :: r, w(8)

    w = this
    r = norm2(w)
  end function

  elemental pure function spin_conjg(this) result(r)
    !! Calculate the complex conjugate of the spin matrix.
    class(spin), intent(in) :: this
    type(spin)              :: r

    r%matrix = conjg(this%matrix)
  end function

  elemental pure function spin_trace(this) result(r)
    !! Calculate the trace of the spin matrix.
    class(spin), intent(in) :: this
    complex(wp)             :: r

    r = this%matrix(1,1) + this%matrix(2,2)
  end function

  pure function spin_inv(this) result(r)
    !! Calculate the inverse of the spin matrix.
    class(spin), intent(in) :: this
    type(spin)              :: r

    r%matrix = matinv2(this%matrix)
  end function

  pure function spin_sum(this) result(r)
    !! Calculate the sum of an array of spin matrices.
    class(spin), intent(in) :: this(:)
    type(spin)              :: r
    integer                 :: n

    do n=1,size(this)
      r % matrix = r % matrix + this(n) % matrix
    end do
  end function
end module
