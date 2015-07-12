! This
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-11
! Updated: 2015-07-11

module module_conductor
  use module_precision
  use module_spin
  use module_state
  implicit none

  ! Class declaration
  type conductor
    type(state), allocatable :: state(:,:) ! Material state as a function of position and energy
    type(state), allocatable :: bcl(:)     ! Boundary condition on the left  as a function of energy
    type(state), allocatable :: bcr(:)     ! Boundary condition on the right as a function of energy
    real(dp),    allocatable :: pos(:)     ! Which positions to consider
    real(dp),    allocatable :: erg(:)     ! Which energies to consider

    contains
    final      :: conductor_destruct      ! Class destructor
  end type

  ! Class constructor
  interface conductor
    module procedure conductor_construct_bcs
  end interface
contains
  pure function conductor_construct_bcs(pos, erg, gap, eps) result(this)
    ! Constructs a conductor object corresponding to a BCS superconductor with a given position and energy range
    type(conductor)                   :: this     ! Conductor object that will be constructed
    real(dp),    intent(in)           :: pos(:)   ! Array of positions where the material is to defined
    real(dp),    intent(in)           :: erg(:)   ! Array of energies where the material is to be defined
    real(dp),    intent(in), optional :: eps      ! Imaginary energy contribution (models inelastic scattering)
    complex(dp), intent(in), optional :: gap      ! Superconducting order parameter

    real(dp)                          :: seps
    complex(dp)                       :: sgap
    integer                           :: n, m

    ! Handle the optional arguments
    if (present(eps)) then
      seps = eps
    else
      seps = 0.001_dp
    end if
   
    if (present(gap)) then
      sgap = gap
    else
      sgap = (1.0_dp,0.0_dp)
    end if

    ! Allocate memory (if necessary)
    if (.not. allocated(this%state)) then
      allocate(this%state(ubound(erg,1), ubound(pos,1)))
      allocate(this%pos(ubound(pos,1)))
      allocate(this%erg(ubound(erg,1)))
      allocate(this%bcl(ubound(erg,1)))
      allocate(this%bcr(ubound(erg,1)))
    end if

    ! Fill the object fields
    this%pos = pos
    this%erg = erg
    forall (n=lbound(erg,1):ubound(erg,1))
      forall (m = lbound(erg,1):ubound(pos,1))
        this%state(n,m) = state( cmplx(erg(n),seps,kind=dp), sgap )
      end forall
      this%bcl(n) = state( cmplx(erg(n),seps,kind=dp), sgap )
      this%bcr(n) = state( cmplx(erg(n),seps,kind=dp), sgap )
    end forall
  end function

  subroutine conductor_destruct(this)
    ! Define the class destructor
    type(conductor) :: this

    ! Deallocate memory (if necessary)
    if(allocated(this%state)) then
      deallocate(this%state)
      deallocate(this%pos)
      deallocate(this%erg)
      deallocate(this%bcl)
      deallocate(this%bcr)
    end if
  end subroutine
end module
