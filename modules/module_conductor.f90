! This module defines the data structure 'conductor', which models the physical state of a metallic 
! conductor as a function of position and energy. The purpose of this structure is mainly to be used
! as a base class for more interesting material classes, such as those that describe superconductors
! and ferromagnets, or to be used in conjunction with such materials in hybrid structures.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-11
! Updated: 2015-07-14

module module_conductor
  use module_precision
  use module_spin
  use module_state
  implicit none

  ! Class declaration
  type conductor
    type(state), allocatable :: state(:,:) ! Physical state as a function of position and energy
    type(state), allocatable :: bcl(:)     ! Boundary condition on the left  as function of energy
    type(state), allocatable :: bcr(:)     ! Boundary condition on the right as function of energy
    real(dp),    allocatable :: pos(:)     ! Discretized position domain 
    real(dp),    allocatable :: erg(:)     ! Discretized energy domain 

    real(dp)  :: interface_left  = inf     ! Interface parameter ζ at the left  interface
    real(dp)  :: interface_right = inf     ! Interface parameter ζ at the right interface
    real(dp)  :: thouless        = 1.0_dp
    real(dp)  :: epsilon         = 0.001_dp

    contains
    final     :: conductor_destruct       ! Class destructor
    procedure :: update => conductor_update
    !TODO class(conductor), pointer :: material_left, material_right
    !TODO procedure Connect_Left / Connect_Right
  end type

  ! Class constructor
  interface conductor
    module procedure conductor_construct_bcs
  end interface
contains
  pure function conductor_construct_bcs(pos, erg, gap, eps) result(this)
    ! Constructs a conductor object corresponding to a BCS superconductor with a given position and energy range
    type(conductor)                   :: this     ! Conductor object that will be constructed
    real(dp),    intent(in)           :: pos(:)   ! Discretized position domain
    real(dp),    intent(in)           :: erg(:)   ! Discretized energy domain
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

  subroutine conductor_update(this)
    ! This subroutine solves the Usadel equation numerically for the given position
    ! and energy range, and updates the current estimate for the state of the system.
    use bvp_m

    class(conductor), intent(inout) :: this
    real(dp)                        :: y(32,size(this%pos))
    type(bvp_sol)                   :: sol
    integer                         :: n, m

    ! Loop over quasiparticle energies
    do n=1,size(this%erg)
      ! Status information
      write (*,'(a,1x,i3,1x,a,1x,i3,1x,a,2x,f8.5)') '[',n,'/',ubound(this%erg,1),']',this%erg(n)

      ! Convert all states at this energy level to real-valued state vectors
      forall (m=1:size(this%pos))
        y(:,m) = this%state(n,m)
      end forall

      ! Initialize the boundary value problem solver (TODO: Provide energy as parameter, compare with LUB)
      sol = bvp_init(32, 16, this%pos, y)

      ! Solve the differential equation
      sol = bvp_solver(sol, ode, bc, method=6, error_control=1, tol=1.0e-6_dp, trace=2)

      ! Use the results to update the current state of the system
      forall (m=1:size(this%pos))
        this%state(n,m) = sol%y(:,m)
      end forall
    end do
  contains
    subroutine ode(t, y, f)
      ! Definition of the differential equation y'=f(t,y)
      real(dp), intent(in)  :: t
      real(dp), intent(in)  :: y(32)
      real(dp), intent(out) :: f(32)

      type(state)           :: s
      type(spin)            :: g, gt, dg, dgt, N, Nt, d2g, d2gt

      ! Convert the state vector to a state object
      s = y

      ! Extract the Riccati parameters
      g   = s%g
      gt  = s%gt
      dg  = s%dg
      dgt = s%dgt

      ! Calculate the normalization matrices
      N   = spin_inv( pauli0 - g*gt )
      Nt  = spin_inv( pauli0 - gt*g )

      ! Calculate the second derivatives of the Riccati parameters according to the Usadel equation TODO fix erg
      d2g  = (-2.0_dp) * dg*Nt*gt*dg - ((0.0_dp,2.0_dp)/this%thouless) * cmplx(this%erg(1),this%epsilon,kind=dp)*g;
      d2gt = (-2.0_dp) * dgt*N*g*dgt - ((0.0_dp,2.0_dp)/this%thouless) * cmplx(this%erg(1),this%epsilon,kind=dp)*gt;
            
      ! Pack the results into a state object
      s%g   = dg
      s%gt  = dgt
      s%dg  = d2g
      s%dgt = d2gt

      ! Convert the state object to a state vector
      f = s
    end subroutine

    subroutine bc(ya,yb,bca,bcb)
      real(dp), intent(in)  :: ya(32)
      real(dp), intent(in)  :: yb(32)
      real(dp), intent(out) :: bca(16)
      real(dp), intent(out) :: bcb(16)

      bca( 1:16) = ya(17:32)
      bcb( 1:16) = yb(17:32)
    end subroutine
  end subroutine
end module
