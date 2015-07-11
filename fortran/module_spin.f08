! This file defines a module containing 
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-10
! Updated: 2015-07-10

module module_spin
  use module_precision
  implicit none

  ! Define a datatype that represents 2×2 complex spin matrices
  type spin
    complex(dp) :: matrix(2,2) = 0
  end type

  interface spin
    module procedure spin_construct_cmatrix, &
                     spin_construct_cvector, &
                     spin_construct_rvector
  end interface
  
  interface operator(*)
    module procedure spin_multl_rscalar, spin_multr_rscalar, &
                     spin_multl_cscalar, spin_multr_cscalar, &
                     spin_multl_cmatrix, spin_multr_cmatrix
  end interface
contains
  function spin_construct_cmatrix(matrix) result(this)
    ! This function constructs a spin object from a 2×2 complex matrix
    type(spin)          :: this
    complex, intent(in) :: matrix(2,2)

    this%matrix = matrix
  end function

  function spin_construct_cvector(vector) result(this)
    ! This function constructs a spin object from a 4×1 complex vector
    type(spin)          :: this
    complex, intent(in) :: vector(4)

    this%matrix = reshape(vector,[2,2])
  end function

  function spin_construct_rvector(vector) result(this)
    ! This function constructs a spin object from a 8×1 real vector
    type(spin)       :: this
    real, intent(in) :: vector(8)

    this%matrix = cmplx( reshape(vector(1:7:2),[2,2]), reshape(vector(2:8:2),[2,2]) )
  end function

  function spin_multl_rscalar(a,b) result(r)
    ! Defines left multiplication of a spin matrix by a real scalar
    type(spin)             :: r
    real,       intent(in) :: a
    type(spin), intent(in) :: b

    r = spin(a * b%matrix)
  end function

  function spin_multr_rscalar(a,b) result(r)
    ! Defines right multiplication of a spin matrix by a real scalar
    type(spin)             :: r
    type(spin), intent(in) :: a
    real,       intent(in) :: b

    r = spin(a%matrix * b)
  end function


  function spin_multl_cscalar(a,b) result(r)
    ! Defines left multiplication of a spin matrix by a complex scalar
    type(spin)             :: r
    complex,    intent(in) :: a
    type(spin), intent(in) :: b

    r = spin(a * b%matrix)
  end function

  function spin_multr_cscalar(a,b) result(r)
    ! Defines right multiplication of a spin matrix by a complex scalar
    type(spin)             :: r
    type(spin), intent(in) :: a
    complex,    intent(in) :: b

    r = spin(a%matrix * b)
  end function

  function spin_multl_cmatrix(a,b) result(r)
    ! Defines left multiplication of a spin matrix by a complex matrix
    type(spin)             :: r
    complex,    intent(in) :: a(2,2)
    type(spin), intent(in) :: b

    r = spin(matmul(a, b%matrix))
  end function

  function spin_multr_cmatrix(a,b) result(r)
    ! Defines right multiplication of a spin matrix by a complex matrix
    type(spin)             :: r
    type(spin), intent(in) :: a
    complex,    intent(in) :: b(2,2)

    r = spin(matmul(a%matrix, b))
  end function
end module
