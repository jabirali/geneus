! This module defines the data type 'spin', which can be used to represent 2×2 complex spin matrices.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-10
! Updated: 2015-07-10

module module_spin
  use module_precision
  implicit none

  ! Class declaration
  type spin
    complex(dp) :: matrix(2,2) = 0
  end type

  ! Class constructor
  interface spin
    module procedure spin_construct_cmatrix, &
                     spin_construct_cvector, &
                     spin_construct_rvector, &
                     spin_construct_spin
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

  ! Exported constants
  type(spin), parameter :: pauli0 = spin(reshape([ ( 1, 0), ( 0, 0), ( 0, 0), ( 1, 0) ], [2,2]))
  type(spin), parameter :: pauli1 = spin(reshape([ ( 0, 0), ( 1, 0), ( 1, 0), ( 0, 0) ], [2,2]))
  type(spin), parameter :: pauli2 = spin(reshape([ ( 0, 0), ( 0,-1), ( 0, 1), ( 0, 0) ], [2,2]))
  type(spin), parameter :: pauli3 = spin(reshape([ ( 1, 0), ( 0, 0), ( 0, 0), (-1, 0) ], [2,2]))
contains
  pure function spin_construct_cmatrix(matrix) result(this)
    ! This function constructs a spin object from a 2×2 complex matrix
    type(spin)          :: this
    complex, intent(in) :: matrix(2,2)

    this%matrix = matrix
  end function

  pure function spin_construct_cvector(vector) result(this)
    ! This function constructs a spin object from a 4×1 complex vector
    type(spin)          :: this
    complex, intent(in) :: vector(4)

    this%matrix = reshape(vector,[2,2])
  end function

  pure function spin_construct_rvector(vector) result(this)
    ! This function constructs a spin object from a 8×1 real vector
    type(spin)       :: this
    real, intent(in) :: vector(8)

    this%matrix = cmplx( reshape(vector(1:7:2),[2,2]), reshape(vector(2:8:2),[2,2]) )
  end function

  pure function spin_construct_spin(other) result(this)
    ! This function constructs a spin object form an existing one
    type(spin)             :: this
    type(spin), intent(in) :: other

    this%matrix = other%matrix
  end function

  pure function spin_multl_rscalar(a,b) result(r)
    ! Defines left multiplication of a spin matrix by a real scalar
    type(spin)             :: r
    real,       intent(in) :: a
    type(spin), intent(in) :: b

    r = spin(a * b%matrix)
  end function

  pure function spin_multr_rscalar(a,b) result(r)
    ! Defines right multiplication of a spin matrix by a real scalar
    type(spin)             :: r
    type(spin), intent(in) :: a
    real,       intent(in) :: b

    r = spin(a%matrix * b)
  end function

  pure function spin_multl_cscalar(a,b) result(r)
    ! Defines left multiplication of a spin matrix by a complex scalar
    type(spin)             :: r
    complex,    intent(in) :: a
    type(spin), intent(in) :: b

    r = spin(a * b%matrix)
  end function

  pure function spin_multr_cscalar(a,b) result(r)
    ! Defines right multiplication of a spin matrix by a complex scalar
    type(spin)             :: r
    type(spin), intent(in) :: a
    complex,    intent(in) :: b

    r = spin(a%matrix * b)
  end function

  pure function spin_multl_cmatrix(a,b) result(r)
    ! Defines left multiplication of a spin matrix by a complex matrix
    type(spin)             :: r
    complex,    intent(in) :: a(2,2)
    type(spin), intent(in) :: b

    r = spin(matmul(a, b%matrix))
  end function

  pure function spin_multr_cmatrix(a,b) result(r)
    ! Defines right multiplication of a spin matrix by a complex matrix
    type(spin)             :: r
    type(spin), intent(in) :: a
    complex,    intent(in) :: b(2,2)

    r = spin(matmul(a%matrix, b))
  end function

  pure function spin_mult_spin(a,b) result(r)
    ! Defines multiplication of two spin matrices
    type(spin)             :: r
    type(spin), intent(in) :: a, b

    r = spin(matmul(a%matrix, b%matrix))
  end function

  pure function spin_addl_rscalar(a,b) result(r)
    ! Defines left addition of a spin matrix and a real scalar
    type(spin)             :: r
    real,       intent(in) :: a
    type(spin), intent(in) :: b

    r = spin(a*pauli0%matrix + b%matrix)
  end function

  pure function spin_addr_rscalar(a,b) result(r)
    ! Defines right addition of a spin matrix and a real scalar
    type(spin)             :: r
    type(spin), intent(in) :: a
    real,       intent(in) :: b

    r = spin(a%matrix + b*pauli0%matrix)
  end function

  pure function spin_addl_cscalar(a,b) result(r)
    ! Defines left addition of a spin matrix and a complex scalar
    type(spin)             :: r
    complex,    intent(in) :: a
    type(spin), intent(in) :: b

    r = spin(a*pauli0%matrix + b%matrix)
  end function

  pure function spin_addr_cscalar(a,b) result(r)
    ! Defines right addition of a spin matrix and a complex scalar
    type(spin)             :: r
    type(spin), intent(in) :: a
    complex,    intent(in) :: b

    r = spin(a%matrix + b*pauli0%matrix)
  end function

  pure function spin_addl_cmatrix(a,b) result(r)
    ! Defines left addition of a spin matrix and a complex matrix
    type(spin)             :: r
    complex,    intent(in) :: a(2,2)
    type(spin), intent(in) :: b

    r = spin(a + b%matrix)
  end function

  pure function spin_addr_cmatrix(a,b) result(r)
    ! Defines right addition of a spin matrix and a complex matrix
    type(spin)             :: r
    type(spin), intent(in) :: a
    complex,    intent(in) :: b(2,2)

    r = spin(a%matrix + b)
  end function

  pure function spin_add_spin(a,b) result(r)
    ! Defines addition of two spin matrices
    type(spin)             :: r
    type(spin), intent(in) :: a, b

    r = spin(a%matrix + b%matrix)
  end function
end module
