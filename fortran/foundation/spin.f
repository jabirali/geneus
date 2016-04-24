!> Author:   Jabir Ali Ouassou
!> Date:     2015-07-10
!> Category: Foundation
!>
!> This module defines the data type 'spin', which represents 2×2 complex matrices in spin space. The module overloads
!> common arithmetic operators to work with the new data type, and defines and exports the Pauli matrices as constants.
!> To make it easier to interact with common differential equation solvers, which often operate on real state vectors,
!> the assignment operator is overloaded in such a way that 'spin' can be easily imported/exported to a real vector(8).
!>
!> @TODO
!>   Replace the remaining methods with generics. In particular, spin_inv and spin_trace should be specific realizations
!>   of a generic function inv and trace from math_m, and spin_print should be replaced by write(formatted) at some point.
!>
!> @TODO
!>   Replace the explicit interface blocks with generic wrappers. The way to do this successfully, seems to be to:
!>    o Define polymorphic procedures, i.e. switch to class(spin),intent(inout) instead of type(spin),intent(in);
!>    o Make these procedures into private type(spin) methods, i.e. import_cscalar => spin_import_cscalar and so on;
!>    o Make public generic interfaces to these procedures, i.e. generic, public :: assignment(=) => import_cscalar, etc.
!>   This approach should work for both the builtin operators, the interfaces conjg etc, and custom ones like trace etc.

module spin_m
  use math_m
  implicit none
  private

  ! Public interface
  public spin
  public assignment(=), operator(+), operator(-), operator(*), operator(/), operator(**), operator(.divl.), operator(.divr.)
  public spin_inv, spin_trace, spin_print, conjg, norm2, sum
  public pauli, pauli0, pauli1, pauli2, pauli3

  ! Type declaration
  type spin
    complex(wp) :: matrix(2,2)      =  0.0_wp        ! Spin matrix
  contains
    procedure   :: inv              => spin_inv      ! Matrix inverse
    procedure   :: trace            => spin_trace    ! Matrix trace
    procedure   :: print            => spin_print    ! Prints the matrix to standard out
  end type

  ! Type constructor
  interface spin
    module procedure spin_construct_rscalar, &
                     spin_construct_cscalar, &
                     spin_construct_cmatrix, &
                     spin_construct_rvector, &
                     spin_construct_spin
  end interface

  ! Assignment operator
  interface assignment(=)
    module procedure spin_import_rscalar, spin_import_cscalar, &
                     spin_import_cmatrix, spin_export_cmatrix, &
                     spin_import_rvector, spin_export_rvector 
  end interface

  ! Exponentiation operator
  interface operator(**)
    module procedure spin_exp_integer
  end interface

  ! Multiplication operator
  interface operator(*)
    module procedure spin_multl_rscalar, spin_multr_rscalar, &
                     spin_multl_cscalar, spin_multr_cscalar, &
                     spin_multl_cmatrix, spin_multr_cmatrix, &
                     spin_mult_spin
  end interface

  ! Division operator
  interface operator(/)
    module procedure spin_div_rscalar, spin_div_cscalar
  end interface

  ! Addition operator
  interface operator(+)
    module procedure spin_addl_rscalar, spin_addr_rscalar, &
                     spin_addl_cscalar, spin_addr_cscalar, &
                     spin_addl_cmatrix, spin_addr_cmatrix, &
                     spin_add_spin
  end interface

  ! Subtraction operator
  interface operator(-)
    module procedure spin_subl_rscalar, spin_subr_rscalar, &
                     spin_subl_cscalar, spin_subr_cscalar, &
                     spin_subl_cmatrix, spin_subr_cmatrix, &
                     spin_sub_spin
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

  ! Matrix division operator (left)
  interface operator(.divl.)
    module procedure spin_divl_spin
  end interface

  ! Matrix division operator (right)
  interface operator(.divr.)
    module procedure spin_divr_spin
  end interface

  ! Exported constants
  type(spin), parameter :: pauli0     = spin(reshape([ ( 1, 0), ( 0, 0), ( 0, 0), ( 1, 0) ], [2,2], order=[2,1]))
  type(spin), parameter :: pauli1     = spin(reshape([ ( 0, 0), ( 1, 0), ( 1, 0), ( 0, 0) ], [2,2], order=[2,1]))
  type(spin), parameter :: pauli2     = spin(reshape([ ( 0, 0), ( 0,-1), ( 0, 1), ( 0, 0) ], [2,2], order=[2,1]))
  type(spin), parameter :: pauli3     = spin(reshape([ ( 1, 0), ( 0, 0), ( 0, 0), (-1, 0) ], [2,2], order=[2,1]))
  type(spin), parameter :: pauli(0:3) = [pauli0, pauli1, pauli2, pauli3]
contains
  pure function spin_construct_rscalar(scalar) result(this)
    ! This function constructs a spin object from a real scalar
    type(spin)             :: this
    real(wp),  intent(in)  :: scalar

    this%matrix = scalar * pauli0%matrix
  end function

  pure function spin_construct_cscalar(scalar) result(this)
    ! This function constructs a spin object from a complex scalar
    type(spin)               :: this
    complex(wp), intent(in)  :: scalar

    this%matrix = scalar * pauli0%matrix
  end function

  pure function spin_construct_cmatrix(matrix) result(this)
    ! This function constructs a spin object from a 2×2 complex matrix
    type(spin)              :: this
    complex(wp), intent(in) :: matrix(2,2)

    this%matrix = matrix
  end function

  pure function spin_construct_rvector(vector) result(this)
    ! This function constructs a spin object from a 8×1 real vector
    type(spin)           :: this
    real(wp), intent(in) :: vector(8)

    this%matrix = cx(reshape(vector(1:7:2),[2,2],order=[2,1]),&
                     reshape(vector(2:8:2),[2,2],order=[2,1]))
  end function

  pure function spin_construct_spin(other) result(this)
    ! This function constructs a spin object form an existing one
    type(spin)             :: this
    type(spin), intent(in) :: other

    this%matrix = other%matrix
  end function

  pure subroutine spin_import_rscalar(this, scalar)
    ! This function assigns a spin object data from a real scalar
    type(spin), intent(out) :: this
    real(wp),   intent(in)  :: scalar

    this%matrix = scalar * pauli0%matrix
  end subroutine

  pure subroutine spin_import_cscalar(this, scalar)
    ! This function assigns a spin object data from a complex scalar
    type(spin),  intent(out) :: this
    complex(wp), intent(in)  :: scalar

    this%matrix = scalar * pauli0%matrix
  end subroutine

  pure subroutine spin_import_cmatrix(this, matrix)
    ! This function assigns a spin object data from a complex matrix
    type(spin),  intent(out) :: this
    complex(wp), intent(in)  :: matrix(2,2)

    this%matrix = matrix
  end subroutine

  pure subroutine spin_import_rvector(this, vector)
    ! This function assigns a spin object data from a real vector
    type(spin), intent(out) :: this
    real(wp),   intent(in)  :: vector(8)

    this%matrix = cx(reshape(vector(1:7:2),[2,2],order=[2,1]),&
                     reshape(vector(2:8:2),[2,2],order=[2,1]))
  end subroutine

  pure subroutine spin_export_cmatrix(matrix, this)
    ! This function assigns a complex matrix from a spin object
    complex(wp), intent(out) :: matrix(2,2)
    type(spin),  intent(in)  :: this

    matrix = this%matrix 
  end subroutine

  pure subroutine spin_export_rvector(vector, this)
    ! This function assigns a real vector from a spin object
    real(wp),  intent(out) :: vector(8)
    type(spin), intent(in) :: this

    vector(1:7:2) = re([ this%matrix(1,:), this%matrix(2,:) ])
    vector(2:8:2) = im([ this%matrix(1,:), this%matrix(2,:) ])
  end subroutine

  pure function spin_exp_integer(a,b) result(r)
    ! Defines exponentiation of a spin matrix by an integer (assumption: positive exponent)
    type(spin)             :: r
    type(spin), intent(in) :: a
    integer,    intent(in) :: b
    integer                :: n

    r = a
    do n=2,b
      r%matrix = r%matrix * a
    end do
  end function

  elemental pure function spin_multl_rscalar(a,b) result(r)
    ! Defines left multiplication of a spin matrix by a real scalar
    type(spin)             :: r
    real(wp),   intent(in) :: a
    type(spin), intent(in) :: b

    r = spin(a * b%matrix)
  end function

  elemental pure function spin_multr_rscalar(a,b) result(r)
    ! Defines right multiplication of a spin matrix by a real scalar
    type(spin)             :: r
    type(spin), intent(in) :: a
    real(wp),   intent(in) :: b

    r = spin(a%matrix * b)
  end function

  elemental pure function spin_multl_cscalar(a,b) result(r)
    ! Defines left multiplication of a spin matrix by a complex scalar
    type(spin)              :: r
    complex(wp), intent(in) :: a
    type(spin),  intent(in) :: b

    r = spin(a * b%matrix)
  end function

  elemental pure function spin_multr_cscalar(a,b) result(r)
    ! Defines right multiplication of a spin matrix by a complex scalar
    type(spin)              :: r
    type(spin),  intent(in) :: a
    complex(wp), intent(in) :: b

    r = spin(a%matrix * b)
  end function

  pure function spin_multl_cmatrix(a,b) result(r)
    ! Defines left multiplication of a spin matrix by a complex matrix
    type(spin)              :: r
    complex(wp), intent(in) :: a(2,2)
    type(spin),  intent(in) :: b
 
    r = spin(matmul(a, b%matrix))
  end function

  pure function spin_multr_cmatrix(a,b) result(r)
    ! Defines right multiplication of a spin matrix by a complex matrix
    type(spin)              :: r
    type(spin),  intent(in) :: a
    complex(wp), intent(in) :: b(2,2)

    r = spin(matmul(a%matrix, b))
  end function

  elemental pure function spin_mult_spin(a,b) result(r)
    ! Defines multiplication of two spin matrices
    type(spin)             :: r
    type(spin), intent(in) :: a, b

    r = spin(matmul(a%matrix, b%matrix))
  end function

  elemental pure function spin_div_rscalar(a,b) result(r)
    ! Defines division of a spin matrix by a real scalar
    type(spin)             :: r
    type(spin), intent(in) :: a
    real(wp),   intent(in) :: b

    r = spin(a%matrix / b)
  end function

  elemental pure function spin_div_cscalar(a,b) result(r)
    ! Defines division of a spin matrix by a real scalar
    type(spin)             :: r
    type(spin), intent(in) :: a
    complex(wp),intent(in) :: b

    r = spin(a%matrix / b)
  end function

  pure function spin_addl_rscalar(a,b) result(r)
    ! Defines left addition of a spin matrix and a real scalar
    type(spin)             :: r
    real(wp),   intent(in) :: a
    type(spin), intent(in) :: b

    r = spin(a*pauli0%matrix + b%matrix)
  end function

  pure function spin_addr_rscalar(a,b) result(r)
    ! Defines right addition of a spin matrix and a real scalar
    type(spin)             :: r
    type(spin), intent(in) :: a
    real(wp),   intent(in) :: b

    r = spin(a%matrix + b*pauli0%matrix)
  end function

  pure function spin_addl_cscalar(a,b) result(r)
    ! Defines left addition of a spin matrix and a complex scalar
    type(spin)              :: r
    complex(wp), intent(in) :: a
    type(spin),  intent(in) :: b

    r = spin(a*pauli0%matrix + b%matrix)
  end function

  pure function spin_addr_cscalar(a,b) result(r)
    ! Defines right addition of a spin matrix and a complex scalar
    type(spin)              :: r
    type(spin),  intent(in) :: a
    complex(wp), intent(in) :: b

    r = spin(a%matrix + b*pauli0%matrix)
  end function

  pure function spin_addl_cmatrix(a,b) result(r)
    ! Defines left addition of a spin matrix and a complex matrix
    type(spin)              :: r
    complex(wp), intent(in) :: a(2,2)
    type(spin),  intent(in) :: b

    r = spin(a + b%matrix)
  end function

  pure function spin_addr_cmatrix(a,b) result(r)
    ! Defines right addition of a spin matrix and a complex matrix
    type(spin)              :: r
    type(spin),  intent(in) :: a
    complex(wp), intent(in) :: b(2,2)

    r = spin(a%matrix + b)
  end function

  pure function spin_add_spin(a,b) result(r)
    ! Defines addition of two spin matrices
    type(spin)             :: r
    type(spin), intent(in) :: a, b

    r = spin(a%matrix + b%matrix)
  end function

  pure function spin_subl_rscalar(a,b) result(r)
    ! Defines left subtraction of a spin matrix and a real scalar
    type(spin)             :: r
    real(wp),   intent(in) :: a
    type(spin), intent(in) :: b

    r = spin(a*pauli0%matrix - b%matrix)
  end function

  pure function spin_subr_rscalar(a,b) result(r)
    ! Defines right subtraction of a spin matrix and a real scalar
    type(spin)             :: r
    type(spin), intent(in) :: a
    real(wp),   intent(in) :: b

    r = spin(a%matrix - b*pauli0%matrix)
  end function

  pure function spin_subl_cscalar(a,b) result(r)
    ! Defines left subtraction of a spin matrix and a complex scalar
    type(spin)              :: r
    complex(wp), intent(in) :: a
    type(spin),  intent(in) :: b

    r = spin(a*pauli0%matrix - b%matrix)
  end function

  pure function spin_subr_cscalar(a,b) result(r)
    ! Defines right subtraction of a spin matrix and a complex scalar
    type(spin)              :: r
    type(spin),  intent(in) :: a
    complex(wp), intent(in) :: b

    r = spin(a%matrix - b*pauli0%matrix)
  end function

  pure function spin_subl_cmatrix(a,b) result(r)
    ! Defines left subtraction of a spin matrix and a complex matrix
    type(spin)              :: r
    complex(wp), intent(in) :: a(2,2)
    type(spin),  intent(in) :: b

    r = spin(a - b%matrix)
  end function

  pure function spin_subr_cmatrix(a,b) result(r)
    ! Defines right subtraction of a spin matrix and a complex matrix
    type(spin)              :: r
    type(spin),  intent(in) :: a
    complex(wp), intent(in) :: b(2,2)

    r = spin(a%matrix - b)
  end function

  pure function spin_sub_spin(a,b) result(r)
    ! Defines subtraction of two spin matrices
    type(spin)             :: r
    type(spin), intent(in) :: a, b

    r = spin(a%matrix - b%matrix)
  end function

  pure function spin_divl_spin(a,b) result(r)
    ! Defines left matrix division of two spin matrices
    type(spin)             :: r
    type(spin), intent(in) :: a, b

    r%matrix = matdivl2(a%matrix, b%matrix)
  end function

  pure function spin_divr_spin(a,b) result(r)
    ! Defines right matrix division of two spin matrices
    type(spin)             :: r
    type(spin), intent(in) :: a, b

    r%matrix = matdivr2(a%matrix, b%matrix)
  end function

  pure function spin_inv(this) result(r)
    ! Calculates the inverse of the spin matrix
    type(spin)              :: r
    class(spin), intent(in) :: this

    r%matrix = matinv2(this%matrix)
  end function

  elemental pure function spin_trace(this) result(r)
    ! Calculates the trace of the spin matrix
    complex(wp)             :: r
    class(spin), intent(in) :: this

    r = this%matrix(1,1) + this%matrix(2,2)
  end function

  elemental pure function spin_norm(this) result(r)
    ! Calculates the Frobenius norm of the spin matrix
    real(wp)                :: r, w(8)
    class(spin), intent(in) :: this

    w = this
    r = norm2(w)
  end function

  pure function spin_sum(this) result(r)
    ! Calculates the sum of an array of spin matrices.
    class(spin), intent(in) :: this(:)
    type(spin)              :: r
    integer                 :: n

    do n=1,size(this)
      r % matrix = r % matrix + this(n) % matrix
    end do
  end function

  elemental pure function spin_conjg(this) result(r)
    ! Calculates the complex conjugate of the spin matrix
    class(spin), intent(in)  :: this
    type(spin)               :: r

    r%matrix = conjg(this%matrix)
  end function

  impure subroutine spin_print(this, title)
    ! Prints the spin matrix to standard out
    class(spin),  intent(in)           :: this 
    character(*), intent(in), optional :: title

    ! Print the name of the matrix if provided
    if(present(title)) then
      print *,' ',title,' = '
    end if

    ! Print the matrix elements
    print '(ss,4x,a,1x,es11.4,1x,a,1x,es11.4,1x,a,5x,es11.4,1x,a,1x,es11.4,1x,a,2x,a)', &
            '⎡',re(this%matrix(1,1)),' +',im(this%matrix(1,1)),'i',                     & 
                re(this%matrix(1,2)),' +',im(this%matrix(1,2)),'i','⎤'
    print '(ss,4x,a,1x,es11.4,1x,a,1x,es11.4,1x,a,5x,es11.4,1x,a,1x,es11.4,1x,a,2x,a)', &
            '⎣',re(this%matrix(2,1)),' +',im(this%matrix(2,1)),'i',                     & 
                re(this%matrix(2,2)),' +',im(this%matrix(2,2)),'i','⎦'

    ! Print extra whitespace
    if(present(title)) then
      print *,''
    end if
  end subroutine
end module
