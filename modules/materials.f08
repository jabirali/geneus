! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-10
! Updated: 2015-07-10

!module mod_spin
!  use mod_prec
!  implicit none
!
!  type Spin
!    complex(dp), dimension(2,2) :: matrix
!    !contains
!    !procedure              :: set   => conductor_set
!    !procedure              :: print => conductor_print
!  end type
!contains
!subroutine cmat2rvec(matrix,vector)
!  ! This function converts a 2×2 complex matrix to an 8×1 real vector. This is used when converting between
!  ! the 2×2 spin matrices that are used in physics, and their numerical representations as real state vectors.
!  complex(dp), dimension(2,2), intent(in)  :: matrix
!  real(dp),    dimension(8),   intent(out) :: vector
!
!  vector = [ real(reshape(matrix,[4])), aimag(reshape(matrix,[4])) ]
!end subroutine
!
!subroutine rvec2cmat(vector,matrix)
!  ! This function converts an 8×1 real vector to a 2×2 complex matrix. This is used when converting between
!  ! the 2×2 spin matrices that are used in physics, and their numerical representations as real state vectors.
!  real(dp),    dimension(8),   intent(in)  :: vector
!  complex(dp), dimension(2,2), intent(out) :: matrix
!
!  matrix = cmplx( reshape(vector(1:4),[2,2]), reshape(vector(5:8),[2,2]) )
!end subroutine 
!end module

program test
  !use materials
  use module_spin
  implicit none

  

  complex(dp), dimension(2,2) :: matrix
  real(dp),    dimension(8)   :: state
  type(spin) :: p 

  matrix = reshape([ (1,2), (3,4), (5,6), (7,8) ], shape(matrix))
  print *,matrix
  state  = matrix
  !call cmat2rvec(state,matrix)
  print *,state
  matrix = state
  !call rvec2cmat(matrix,state)
  print *,matrix

  !p = real([ 1., 2., 3., 4., 5., 6., 7., 8. ], kind=dp)
  call p%set_state([ 1d0, 2d0, 3d0, 4d0, 5d0, 6d0, 7d0, 8d0 ])
  print *,p
  !p%set_state([ 1., 2., 3., 4., 5., 6., 7., 8. ])
  print *,pauli(3,:,:)

end program 



  !type(Conductor) :: c
  !call c%set([ 1., 2., 3., 4., 5. ])
  !call c%print


!module materials
!  use mod_prec
!  implicit none
!
!  type Conductor
!    complex(dp)            :: test
!    real(dp), dimension(5) :: content
!    contains
!    procedure              :: set   => conductor_set
!    procedure              :: print => conductor_print
!  end type Conductor
!
!  contains
!    subroutine conductor_set(this, content)
!      class(Conductor),   intent(out) :: this
!      real, dimension(:), intent(in)  :: content
!
!      this%test = complex(1,3)
!
!      this%content = content
!    end subroutine conductor_set
!
!    subroutine conductor_print(this)
!      class(Conductor), intent(in)  :: this
!
!      print *,[ real( [ complex(1,2), complex(3,4) ] ), aimag( [ complex(1,2), complex(3,4) ] ) ]
!
!      print *,'Content is: ',this%content
!    end subroutine conductor_print
!end module materials

