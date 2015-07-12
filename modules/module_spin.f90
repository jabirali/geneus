! This module defines the data type 'spin', which can be used to represent 2×2 complex spin matrices,
! and overloads some common arithmetic operators to work for the new data type. The module then uses
! the data type to define and export the Pauli matrices ('pauli0', 'pauli1', 'pauli2', and 'pauli3').
! To make it easier to interact with common ODE solvers, which often operate on real state vectors,
! the assignment operator is also overloaded in such a way that 'spin' becomes isomorphic to real(8).
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-10
! Updated: 2015-07-11

module module_spin
  use module_precision
  implicit none

  ! Class declaration
  type spin
    complex(dp) :: matrix(2,2) = 0                   ! Stores the spin matrix
    contains
    procedure   :: trace            => spin_trace    ! Calculates the trace of the spin matrix
    procedure   :: inv              => spin_inv      ! Calculates the inverse of the spin matrix
    procedure   :: print            => spin_print    ! Prints the spin matrix to standard out
   !generic     :: write(formatted) => spin_write    ! Modifies the output format (used by 'write' and 'print')
   !generic     :: read(formatted)  => spin_write    ! Modifies the input format (used by 'read')
  end type

  ! Class constructor
  interface spin
    module procedure spin_construct_cmatrix, &
                     spin_construct_cvector, &
                     spin_construct_rvector, &
                     spin_construct_spin
  end interface

  ! Assignment operator
  interface assignment(=)
    module procedure spin_import_rscalar, spin_import_cscalar, &
                     spin_import_cmatrix, spin_export_cmatrix, &
                     spin_import_cvector, spin_export_cvector, &
                     spin_import_rvector, spin_export_rvector 
  end interface
  
  ! Multiplication operator
  interface operator(*)
    module procedure spin_multl_rscalar, spin_multr_rscalar, &
                     spin_multl_cscalar, spin_multr_cscalar, &
                     spin_multl_cmatrix, spin_multr_cmatrix, &
                     spin_mult_spin
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

  ! Matrix division operator (left)
  interface operator(.divl.)
    module procedure spin_divl_spin
  end interface

  ! Matrix division operator (right)
  interface operator(.divr.)
    module procedure spin_divr_spin
  end interface

  ! Exported constants
  type(spin), parameter :: pauli0 = spin(reshape([ ( 1, 0), ( 0, 0), ( 0, 0), ( 1, 0) ], [2,2]))
  type(spin), parameter :: pauli1 = spin(reshape([ ( 0, 0), ( 1, 0), ( 1, 0), ( 0, 0) ], [2,2]))
  type(spin), parameter :: pauli2 = spin(reshape([ ( 0, 0), ( 0, 1), ( 0,-1), ( 0, 0) ], [2,2]))
  type(spin), parameter :: pauli3 = spin(reshape([ ( 1, 0), ( 0, 0), ( 0, 0), (-1, 0) ], [2,2]))
contains
  pure function spin_construct_cmatrix(matrix) result(this)
    ! This function constructs a spin object from a 2×2 complex matrix
    type(spin)              :: this
    complex(dp), intent(in) :: matrix(2,2)

    this%matrix = matrix
  end function

  pure function spin_construct_cvector(vector) result(this)
    ! This function constructs a spin object from a 4×1 complex vector
    type(spin)              :: this
    complex(dp), intent(in) :: vector(4)

    this%matrix = reshape(vector,[2,2])
  end function

  pure function spin_construct_rvector(vector) result(this)
    ! This function constructs a spin object from a 8×1 real vector
    type(spin)           :: this
    real(dp), intent(in) :: vector(8)

    this%matrix = cmplx( reshape(vector(1:7:2),[2,2]), reshape(vector(2:8:2),[2,2]) )
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
    real(dp),   intent(in)  :: scalar

    this%matrix = scalar * pauli0%matrix
  end subroutine

  pure subroutine spin_import_cscalar(this, scalar)
    ! This function assigns a spin object data from a complex scalar
    type(spin),  intent(out) :: this
    complex(dp), intent(in)  :: scalar

    this%matrix = scalar * pauli0%matrix
  end subroutine

  pure subroutine spin_import_cmatrix(this, matrix)
    ! This function assigns a spin object data from a complex matrix
    type(spin),  intent(out) :: this
    complex(dp), intent(in)  :: matrix(2,2)

    this%matrix = matrix
  end subroutine

  pure subroutine spin_import_cvector(this, vector)
    ! This function assigns a spin object data from a complex vector
    type(spin),  intent(out) :: this
    complex(dp), intent(in)  :: vector(4)

    this%matrix = reshape(vector,[2,2])
  end subroutine

  pure subroutine spin_import_rvector(this, vector)
    ! This function assigns a spin object data from a real vector
    type(spin), intent(out) :: this
    real(dp),   intent(in)  :: vector(8)

    this%matrix = cmplx( reshape(vector(1:7:2),[2,2]), reshape(vector(2:8:2),[2,2]) )
  end subroutine

  pure subroutine spin_export_cmatrix(matrix, this)
    ! This function assigns a complex matrix from a spin object
    complex(dp), intent(out) :: matrix(2,2)
    type(spin),  intent(in)  :: this

    matrix = this%matrix 
  end subroutine

  pure subroutine spin_export_cvector(vector, this)
    ! This function assigns a complex vector from a spin object
    complex(dp), intent(out) :: vector(4)
    type(spin),  intent(in)  :: this

    vector = reshape(this%matrix,[4])
  end subroutine

  pure subroutine spin_export_rvector(vector, this)
    ! This function assigns a real vector from a spin object
    real(dp),  intent(out) :: vector(8)
    type(spin), intent(in) :: this

    vector(1:7:2) =  real(reshape(this%matrix,[4]))
    vector(2:8:2) = aimag(reshape(this%matrix,[4]))
  end subroutine

  pure function spin_multl_rscalar(a,b) result(r)
    ! Defines left multiplication of a spin matrix by a real scalar
    type(spin)             :: r
    real(dp),   intent(in) :: a
    type(spin), intent(in) :: b

    r = spin(a * b%matrix)
  end function

  pure function spin_multr_rscalar(a,b) result(r)
    ! Defines right multiplication of a spin matrix by a real scalar
    type(spin)             :: r
    type(spin), intent(in) :: a
    real(dp),   intent(in) :: b

    r = spin(a%matrix * b)
  end function

  pure function spin_multl_cscalar(a,b) result(r)
    ! Defines left multiplication of a spin matrix by a complex scalar
    type(spin)             :: r
    complex(dp), intent(in) :: a
    type(spin),  intent(in) :: b

    r = spin(a * b%matrix)
  end function

  pure function spin_multr_cscalar(a,b) result(r)
    ! Defines right multiplication of a spin matrix by a complex scalar
    type(spin)              :: r
    type(spin),  intent(in) :: a
    complex(dp), intent(in) :: b

    r = spin(a%matrix * b)
  end function

  pure function spin_multl_cmatrix(a,b) result(r)
    ! Defines left multiplication of a spin matrix by a complex matrix
    type(spin)              :: r
    complex(dp), intent(in) :: a(2,2)
    type(spin),  intent(in) :: b
 
    ! TODO: Compare this to 'zgemm'
    r = spin(matmul(a, b%matrix))
  end function

  pure function spin_multr_cmatrix(a,b) result(r)
    ! Defines right multiplication of a spin matrix by a complex matrix
    type(spin)              :: r
    type(spin),  intent(in) :: a
    complex(dp), intent(in) :: b(2,2)

    ! TODO: Compare this to 'zgemm'
    r = spin(matmul(a%matrix, b))
  end function

  pure function spin_mult_spin(a,b) result(r)
    ! Defines multiplication of two spin matrices
    type(spin)             :: r
    type(spin), intent(in) :: a, b

    ! TODO: Compare this to 'zgemm'
    r = spin(matmul(a%matrix, b%matrix))
  end function

  pure function spin_addl_rscalar(a,b) result(r)
    ! Defines left addition of a spin matrix and a real scalar
    type(spin)             :: r
    real(dp),   intent(in) :: a
    type(spin), intent(in) :: b

    r = spin(a*pauli0%matrix + b%matrix)
  end function

  pure function spin_addr_rscalar(a,b) result(r)
    ! Defines right addition of a spin matrix and a real scalar
    type(spin)             :: r
    type(spin), intent(in) :: a
    real(dp),   intent(in) :: b

    r = spin(a%matrix + b*pauli0%matrix)
  end function

  pure function spin_addl_cscalar(a,b) result(r)
    ! Defines left addition of a spin matrix and a complex scalar
    type(spin)              :: r
    complex(dp), intent(in) :: a
    type(spin),  intent(in) :: b

    r = spin(a*pauli0%matrix + b%matrix)
  end function

  pure function spin_addr_cscalar(a,b) result(r)
    ! Defines right addition of a spin matrix and a complex scalar
    type(spin)              :: r
    type(spin),  intent(in) :: a
    complex(dp), intent(in) :: b

    r = spin(a%matrix + b*pauli0%matrix)
  end function

  pure function spin_addl_cmatrix(a,b) result(r)
    ! Defines left addition of a spin matrix and a complex matrix
    type(spin)              :: r
    complex(dp), intent(in) :: a(2,2)
    type(spin),  intent(in) :: b

    r = spin(a + b%matrix)
  end function

  pure function spin_addr_cmatrix(a,b) result(r)
    ! Defines right addition of a spin matrix and a complex matrix
    type(spin)              :: r
    type(spin),  intent(in) :: a
    complex(dp), intent(in) :: b(2,2)

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
    real(dp),   intent(in) :: a
    type(spin), intent(in) :: b

    r = spin(a*pauli0%matrix - b%matrix)
  end function

  pure function spin_subr_rscalar(a,b) result(r)
    ! Defines right subtraction of a spin matrix and a real scalar
    type(spin)             :: r
    type(spin), intent(in) :: a
    real(dp),   intent(in) :: b

    r = spin(a%matrix - b*pauli0%matrix)
  end function

  pure function spin_subl_cscalar(a,b) result(r)
    ! Defines left subtraction of a spin matrix and a complex scalar
    type(spin)              :: r
    complex(dp), intent(in) :: a
    type(spin),  intent(in) :: b

    r = spin(a*pauli0%matrix - b%matrix)
  end function

  pure function spin_subr_cscalar(a,b) result(r)
    ! Defines right subtraction of a spin matrix and a complex scalar
    type(spin)              :: r
    type(spin),  intent(in) :: a
    complex(dp), intent(in) :: b

    r = spin(a%matrix - b*pauli0%matrix)
  end function

  pure function spin_subl_cmatrix(a,b) result(r)
    ! Defines left subtraction of a spin matrix and a complex matrix
    type(spin)              :: r
    complex(dp), intent(in) :: a(2,2)
    type(spin),  intent(in) :: b

    r = spin(a - b%matrix)
  end function

  pure function spin_subr_cmatrix(a,b) result(r)
    ! Defines right subtraction of a spin matrix and a complex matrix
    type(spin)              :: r
    type(spin),  intent(in) :: a
    complex(dp), intent(in) :: b(2,2)

    r = spin(a%matrix - b)
  end function

  pure function spin_sub_spin(a,b) result(r)
    ! Defines subtraction of two spin matrices
    type(spin)             :: r
    type(spin), intent(in) :: a, b

    r = spin(a%matrix - b%matrix)
  end function

  function spin_divl_spin(a,b) result(r)
    ! Defines left matrix division of two spin matrices
    type(spin)             :: r
    type(spin), intent(in) :: a, b
    
    integer, parameter     :: n = 2      ! Size of matrices A and B
    complex(dp)            :: la(n,n)    ! Copy of A sent to LAPACK
    complex(dp)            :: lb(n,n)    ! Copy of B sent to LAPACK
    integer                :: ipiv(n)    ! Pivot indices 
    integer                :: info       ! Nonzero if any errors occured

    ! Copy the contents of A and B to mutable matrices
    la = a
    lb = b

    ! Call LAPACK to solve the equation AX=B
    call zgesv( n, n, la, n, ipiv, lb, n, info )

    ! Return the result as a spin object
    r = lb
  end function

  function spin_divr_spin(a,b) result(r)
    ! Defines right matrix division of two spin matrices
    type(spin)             :: r
    type(spin), intent(in) :: a, b

    r = a * b%inv()
  end function


  pure function spin_trace(this) result(r)
    ! Calculates the trace of the spin matrix
    complex(dp)             :: r
    class(spin), intent(in) :: this

    r = this%matrix(1,1) + this%matrix(2,2)
  end function

  function spin_inv(this) result(r)
    ! Calculates the inverse of the spin matrix
    type(spin)              :: r
    class(spin), intent(in) :: this

    integer, parameter      :: n = 2    ! Size of the matrix
    integer, parameter      :: m = 64*n ! Size of the work array
    complex(dp)             :: a(n,n)   ! Local copy of the matrix
    complex(dp)             :: work(m)  ! Local work array
    integer                 :: ipiv(n)  ! Pivot indices
    integer                 :: info     ! Nonzero if any errors occured

    ! Copy the contents of this to a mutable matrix
    a = this

    ! Call LAPACK to perform an LU factorization of the matrix
    call zgetrf( n, n, a, n, ipiv, info )

    ! Call LAPACK to invert the matrix
    call zgetri( n, a, n, ipiv, work, m, info )

    ! Return the result as a spin object
    r = a
  end function

  subroutine spin_print(this, title)
    ! Prints the spin matrix to stdout
    class(spin),  intent(in)           :: this 
    character(*), intent(in), optional :: title

    ! Print the name of the matrix if provided
    if(present(title)) then
      print *,title,' = '
    end if

    ! Print the matrix elements
    write(*,'(ss,4x,a,2x,es11.4,1x,a,1x,es11.4,1x,a,4x,es11.4,1x,a,1x,es11.4,1x,a,2x,a)') &
            '⎛',real(this%matrix(1,1)),'+',aimag(this%matrix(1,1)),'i',                & 
                real(this%matrix(1,2)),'+',aimag(this%matrix(1,2)),'i','⎞'
    write(*,'(ss,4x,a,2x,es11.4,1x,a,1x,es11.4,1x,a,4x,es11.4,1x,a,1x,es11.4,1x,a,2x,a)') &
            '⎝',real(this%matrix(2,1)),'+',aimag(this%matrix(2,1)),'i',                & 
                real(this%matrix(2,2)),'+',aimag(this%matrix(2,2)),'i','⎠'
  end subroutine

 !subroutine spin_read(this, unit, iotype, v_list, iostat, iomsg)
 !  ! This method is used to overload the builtin 'read'
 !  class(spin),  intent(inout) :: this
 !  integer,      intent(in   ) :: unit
 !  character(*), intent(in   ) :: iotype
 !  integer,      intent(in   ) :: v_list(:)
 !  integer,      intent(  out) :: iostat
 !  character(*), intent(inout) :: iomsg
 !
 !  read(unit, fmt=*, iostat=iostat, iomsg=iomsg) this%matrix
 !end subroutine

 !subroutine spin_write(this, unit, iotype, v_list, iostat, iomsg)
 !  ! This method is used to overload the builtin 'write'
 !  class(spin),  intent(in   ) :: this
 !  integer,      intent(in   ) :: unit
 !  character(*), intent(in   ) :: iotype
 !  integer,      intent(  out) :: v_list(:)
 !  integer,      intent(  out) :: iostat
 !  character(*), intent(inout) :: iomsg
 !
 !  write(unit, fmt=*, iostat=iostat, iomsg=iomsg) this%matrix
 !end subroutine
end module
