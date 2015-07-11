! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-10
! Updated: 2015-07-10

module module_spin
  use module_precision
  implicit none

  complex(dp), dimension(3,2,2), parameter ::         &
    pauli = reshape([ ( 0, 0),  ( 0, 0), ( 1, 0),     &
                      ( 1, 0),  ( 0,-1), ( 0, 0),     &
                      ( 1, 0),  ( 0, 1), ( 0, 0),     &
                      ( 0, 0),  ( 0, 0), (-1, 0)  ], [3,2,2])

  complex(dp), dimension(2,2), parameter :: paulix = pauli(1,:,:)
  complex(dp), dimension(2,2), parameter :: pauliy = pauli(2,:,:)
  complex(dp), dimension(2,2), parameter :: pauliz = pauli(3,:,:)

  type spin
    complex(dp), dimension(2,2) :: matrix
  contains
    !generic   :: assignment(=) => set_state
    procedure :: set_state     => spin_set_state
    procedure :: get_state     => spin_get_state
  end type

  interface assignment(=)
    module procedure cmat2rvec, rvec2cmat
  end interface

  !interface operator(*)
  !  module procedure multest
  !end interface

contains
  subroutine ass(a,b)
    type(spin), intent(out) :: a
    real(dp),   intent(in)  :: b(8)

    call spin_set_state(a,b)
  end subroutine
  subroutine spin_set_state(this,state)
    ! This method initializes the complex(2,2) spin matrix from a real(8) state vector
    class(spin), intent(out) :: this
    real(dp),   intent(in)   :: state(8)

    this%matrix = cmplx( reshape(state(1:4),[2,2]), reshape(state(5:8),[2,2]) )
  end subroutine

  function multest(a,b) result(c)
    implicit none
    complex(dp), dimension(3,2,2), intent(in)  :: a,b
    complex(dp), dimension(2,2)                :: c

    c = matmul(a(1,:,:), b(1,:,:)) + matmul(a(2,:,:), b(2,:,:)) + matmul(a(3,:,:), b(3,:,:))
  end function

  subroutine spin_get_state(this,state)
    ! This method converts the complex(2,2) spin matrix to a real(8) state vector
    class(spin), intent(in) :: this
    real(dp),   intent(out) :: state(8)

    state = [ real(reshape(this%matrix,[4])), aimag(reshape(this%matrix,[4])) ]
  end subroutine

  subroutine cmat2rvec(vector,matrix)
    ! This function converts a 2×2 complex matrix to an 8×1 real vector. This is used when converting between
    ! the 2×2 spin matrices that are used in physics, and their numerical representations as real state vectors.
    complex(dp), dimension(2,2), intent(in)  :: matrix
    real(dp),    dimension(8),   intent(out) :: vector
  
    vector = [ real(reshape(matrix,[4])), aimag(reshape(matrix,[4])) ]
  end subroutine
  
  subroutine rvec2cmat(matrix,vector)
    ! This function converts an 8×1 real vector to a 2×2 complex matrix. This is used when converting between
    ! the 2×2 spin matrices that are used in physics, and their numerical representations as real state vectors.
    real(dp),    dimension(8),   intent(in)  :: vector
    complex(dp), dimension(2,2), intent(out) :: matrix
  
    matrix = cmplx( reshape(vector(1:4),[2,2]), reshape(vector(5:8),[2,2]) )
  end subroutine 
end module

