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
    ! User-defined parameters: these control the physical characteristics of the material
    real(dp)  :: thouless        = 1.000_dp        ! Thouless energy of the material
    real(dp)  :: scattering      = 0.001_dp        ! Imaginary energy contribution due to inelastic scattering
    real(dp)  :: conductance_a   = 1.0_dp/3.0_dp   ! Tunneling conductance of the left interface  (relative to the bulk conductance of this material)
    real(dp)  :: conductance_b   = 0.0_dp          ! Tunneling conductance of the right interface (relative to the bulk conductance of this material)

    ! Core internal structure: these essential variables store the physical state of the material
    type(state), allocatable :: state(:,:)         ! Physical state as a function of position and energy
    type(state), allocatable :: state_a(:)         ! Boundary condition on the left  as function of energy
    type(state), allocatable :: state_b(:)         ! Boundary condition on the right as function of energy
    real(dp),    allocatable :: pos(:)             ! Discretized position domain 
    real(dp),    allocatable :: erg(:)             ! Discretized energy domain 

    ! Temp internal structure: these private variables are only used by internal subroutines
    type(state), private     :: work_state_a
    type(state), private     :: work_state_b
    complex(dp), private     :: work_erg

    contains
    final     :: conductor_destruct       ! Class destructor
    procedure :: update => conductor_update
    !TODO class(conductor), pointer :: material_a, material_right
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
      allocate(this%state(size(erg), size(pos)))
      allocate(this%pos(size(pos)))
      allocate(this%erg(size(erg)))
      allocate(this%state_a(size(erg)))
      allocate(this%state_b(size(erg)))
    end if

    ! Fill the object fields
    this%pos = pos
    this%erg = erg
    forall (n=lbound(erg,1):ubound(erg,1))
      forall (m = lbound(erg,1):ubound(pos,1))
        this%state(n,m) = state( cmplx(erg(n),seps,kind=dp), sgap )
      end forall
      this%state_a(n) = state( cmplx(erg(n),seps,kind=dp), sgap )
      this%state_b(n) = state( cmplx(erg(n),seps,kind=dp), sgap )
    end forall
  end function

  pure subroutine conductor_destruct(this)
    ! Define the class destructor
    type(conductor), intent(inout) :: this

    ! Deallocate memory (if necessary)
    if(allocated(this%state)) then
      deallocate(this%state)
      deallocate(this%state_a)
      deallocate(this%state_b)
      deallocate(this%pos)
      deallocate(this%erg)
    end if
  end subroutine

  subroutine conductor_update(this)
    ! This subroutine solves the Usadel equation numerically for the given position
    ! and energy range, and updates the current estimate for the state of the system.
    use bvp_m

    class(conductor), intent(inout) :: this
    real(dp)                        :: u(32,size(this%pos))
    type(bvp_sol)                   :: sol
    integer                         :: n, m

    ! Loop over quasiparticle energies
    do n=1,size(this%erg)
      ! Status information
      write (*,'(a,1x,i3,1x,a,1x,i3,1x,a,2x,f8.5)') '[',n,'/',size(this%erg),']',this%erg(n)

      ! Convert all states at this energy level to real-valued state vectors
      forall (m=1:size(this%pos))
        u(:,m) = this%state(n,m)
      end forall

      ! Copy the boundary conditions and energy to temporary work variables
      this%work_state_a = this%state_a(n)
      this%work_state_b = this%state_b(n)
      this%work_erg = cmplx(this%erg(n)/this%thouless, this%scattering/this%thouless, kind=dp)

      ! Initialize the boundary value problem solver
      sol = bvp_init(32, 16, this%pos, u)

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
    subroutine ode(z, u, f)
      ! Definition of the differential equation u'=f(z,u)
      real(dp), intent(in)  :: z
      real(dp), intent(in)  :: u(32)
      real(dp), intent(out) :: f(32)

      type(spin)            :: g, gt, dg, dgt, d2g, d2gt, N, Nt

      ! Extract the Riccati parameters
      g   = u( 1: 8)
      gt  = u( 9:16)
      dg  = u(17:24)
      dgt = u(25:32)

      ! Calculate the normalization matrices
      N   = spin_inv( pauli0 - g*gt )
      Nt  = spin_inv( pauli0 - gt*g )

      ! Calculate the second derivatives of the Riccati parameters according to the Usadel equation
      d2g  = (-2.0_dp)*dg*Nt*gt*dg - (0.0_dp,2.0_dp)*this%work_erg*g
      d2gt = (-2.0_dp)*dgt*N*g*dgt - (0.0_dp,2.0_dp)*this%work_erg*gt
            
      ! Pack the results into a state object
      f( 1: 8) = dg
      f( 9:16) = dgt
      f(17:24) = d2g
      f(25:32) = d2gt
    end subroutine

    subroutine bc(ua, ub, bca, bcb)
      real(dp), intent(in)  :: ua(32)
      real(dp), intent(in)  :: ub(32)
      real(dp), intent(out) :: bca(16)
      real(dp), intent(out) :: bcb(16)

      type(spin)            :: g0, gt0, dg0, dgt0, N0, Nt0
      type(spin)            :: g1, gt1, dg1, dgt1
      type(spin)            :: g2, gt2, dg2, dgt2
      type(spin)            :: g3, gt3, dg3, dgt3, N3, Nt3

      ! State at the left end of the material
      g1   = ua( 1: 8)
      gt1  = ua( 9:16)
      dg1  = ua(17:24)
      dgt1 = ua(25:32)

      ! State at the right end of the material
      g2   = ub( 1: 8)
      gt2  = ub( 9:16)
      dg2  = ub(17:24)
      dgt2 = ub(25:32)

      ! Left boundary condition
      if (this%conductance_a > 0.0_dp) then
        ! State in the material to the left
        g0   = this%work_state_a%g
        gt0  = this%work_state_a%gt
        dg0  = this%work_state_a%dg
        dgt0 = this%work_state_a%dgt

        ! Calculate the normalization matrices
        N0  = spin_inv( pauli0 - g0*gt0 )
        Nt0 = spin_inv( pauli0 - gt0*g0 )

        ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
        bca(1: 8) = dg1  - this%conductance_a*( pauli0 - g1*gt0 )*N0*(  g1  - g0  )
        bca(9:16) = dgt1 - this%conductance_a*( pauli0 - gt1*g0 )*Nt0*( gt1 - gt0 )
      else
        ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
        bca(1: 8) = dg1
        bca(9:16) = dgt1
      end if
        
      ! Right boundary condition
      if (this%conductance_b > 0.0_dp) then
        ! State in the material to the right
        g3   = this%work_state_b%g;
        gt3  = this%work_state_b%gt;
        dg3  = this%work_state_b%dg;
        dgt3 = this%work_state_b%dgt;

        ! Calculate the normalization matrices
        N3  = spin_inv( pauli0 - g3*gt3 );
        Nt3 = spin_inv( pauli0 - gt3*g3 );
        
        ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
        bcb(1: 8) = dg2  - this%conductance_b*( pauli0 - g2*gt3 )*N3*(  g3  - g2  )
        bcb(9:16) = dgt2 - this%conductance_b*( pauli0 - gt2*g3 )*Nt3*( gt3 - gt2 )
      else
        ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
        bcb(1: 8) = dg2
        bcb(9:16) = dgt2
      end if
    end subroutine
  end subroutine
end module
