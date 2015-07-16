! This module defines the data structure 'conductor', which models the physical state of a metallic 
! conductor as a function of position and energy. The purpose of this structure is mainly to be used
! as a base class for more interesting material classes, such as those that describe superconductors
! and ferromagnets, or to be used in conjunction with such materials to describe hybrid structures.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-11
! Updated: 2015-07-15

module module_conductor
  use module_precision
  use module_spin
  use module_state
  implicit none

  ! Class declaration
  type conductor
    type(state), allocatable :: state(:,:) ! Physical state as a function of position and energy
    type(state), allocatable :: state_left(:)     ! Boundary condition on the left  as function of energy
    type(state), allocatable :: state_right(:)     ! Boundary condition on the right as function of energy
    real(dp),    allocatable :: pos(:)     ! Discretized position domain 
    real(dp),    allocatable :: erg(:)     ! Discretized energy domain 

    real(dp)  :: interface_left  = 3.0_dp     ! Interface parameter ζ at the left  interface
    real(dp)  :: interface_right = inf     ! Interface parameter ζ at the right interface
    real(dp)  :: thouless        = 1.000_dp
    real(dp)  :: epsilon         = 0.001_dp

    !real(dp),    private, allocatable :: work_state(:,:)
    !real(dp),    private :: work_state(:,:)
    type(state), private :: work_state_left
    type(state), private :: work_state_right
    complex(dp), private :: work_energy

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
      allocate(this%state_left(ubound(erg,1)))
      allocate(this%state_right(ubound(erg,1)))
    end if

    ! Fill the object fields
    this%pos = pos
    this%erg = erg
    forall (n=lbound(erg,1):ubound(erg,1))
      forall (m = lbound(erg,1):ubound(pos,1))
        this%state(n,m) = state( cmplx(erg(n),seps,kind=dp), sgap )
      end forall
      this%state_left(n) = state( cmplx(erg(n),seps,kind=dp), sgap )
      this%state_right(n) = state( cmplx(erg(n),seps,kind=dp), sgap )
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
      deallocate(this%state_left)
      deallocate(this%state_right)
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
      this%work_state_left = this%state_left(n)
      this%work_state_right = this%state_right(n)

      ! UNKNOWN parameters in p. Replace with global variables in class.
      this%work_energy = cmplx(this%erg(n)/this%thouless, this%epsilon/this%thouless, kind=dp)

      ! Initialize the boundary value problem solver (TODO: Provide energy as parameter, compare with LUB)
      sol = bvp_init(32, 16, this%pos, y)

      ! Solve the differential equation
      sol = bvp_solver(sol, ode, bc, method=6, error_control=1, tol=1.0e-6_dp, trace=0)

      ! Use the results to update the current state of the system
      forall (m=1:size(this%pos))
        this%state(n,m) = sol%y(:,m)
      end forall

    end do
    ! Clean up
    call bvp_terminate(sol)
  contains
    subroutine ode(t, y, f)
      ! Definition of the differential equation y'=f(t,y)
      real(dp), intent(in)  :: t
      real(dp), intent(in)  :: y(32)
      real(dp), intent(out) :: f(32)

      type(spin)            :: g, gt, dg, dgt, d2g, d2gt, N, Nt

      ! Extract the Riccati parameters
      g   = y( 1: 8)
      gt  = y( 9:16)
      dg  = y(17:24)
      dgt = y(25:32)

      ! Calculate the normalization matrices
      N   = spin_inv( pauli0 - g*gt )
      Nt  = spin_inv( pauli0 - gt*g )

      ! Calculate the second derivatives of the Riccati parameters according to the Usadel equation
      d2g  = (-2.0_dp)*dg*Nt*gt*dg - (0.0_dp,2.0_dp)*this%work_energy*g
      d2gt = (-2.0_dp)*dgt*N*g*dgt - (0.0_dp,2.0_dp)*this%work_energy*gt
            
      ! Pack the results into a state object
      f( 1: 8) = dg
      f( 9:16) = dgt
      f(17:24) = d2g
      f(25:32) = d2gt
    end subroutine

    subroutine bc(ya,yb,bca,bcb)
      real(dp), intent(in)  :: ya(32)
      real(dp), intent(in)  :: yb(32)
      real(dp), intent(out) :: bca(16)
      real(dp), intent(out) :: bcb(16)

      type(spin)            :: g0, gt0, dg0, dgt0, N0, Nt0
      type(spin)            :: g1, gt1, dg1, dgt1
      type(spin)            :: g2, gt2, dg2, dgt2
      type(spin)            :: g3, gt3, dg3, dgt3, N3, Nt3

      ! State at the left end of the material
      g1   = ya( 1: 8)
      gt1  = ya( 9:16)
      dg1  = ya(17:24)
      dgt1 = ya(25:32)

      ! State at the right end of the material
      g2   = yb( 1: 8)
      gt2  = yb( 9:16)
      dg2  = yb(17:24)
      dgt2 = yb(25:32)

      ! Left boundary condition
      if (this%interface_left < inf) then
        ! State in the material to the left
        g0   = this%work_state_left%g
        gt0  = this%work_state_left%gt
        dg0  = this%work_state_left%dg
        dgt0 = this%work_state_left%dgt

        ! Calculate the normalization matrices
        N0  = spin_inv( pauli0 - g0*gt0 )
        Nt0 = spin_inv( pauli0 - gt0*g0 )

        ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
        bca(1: 8) = dg1  - ( pauli0 - g1*gt0 )*N0*(  g1  - g0  )/this%interface_left
        bca(9:16) = dgt1 - ( pauli0 - gt1*g0 )*Nt0*( gt1 - gt0 )/this%interface_left
      else
        ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
        bca(1: 8) = dg1
        bca(9:16) = dgt1
      end if
        
      ! Right boundary condition
      if (this%interface_right < inf) then
        ! State in the material to the right
        g3   = this%work_state_right%g;
        gt3  = this%work_state_right%gt;
        dg3  = this%work_state_right%dg;
        dgt3 = this%work_state_right%dgt;

        ! Calculate the normalization matrices
        N3  = spin_inv( pauli0 - g3*gt3 );
        Nt3 = spin_inv( pauli0 - gt3*g3 );
        
        ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
        bcb(1: 8) = dg2  - ( pauli0 - g2*gt3 )*N3*(  g3  - g2  )/this%interface_right
        bcb(9:16) = dgt2 - ( pauli0 - gt2*g3 )*Nt3*( gt3 - gt2 )/this%interface_right
      else
        ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
        bcb(1: 8) = dg2
        bcb(9:16) = dgt2
      end if
    end subroutine
  end subroutine
end module
