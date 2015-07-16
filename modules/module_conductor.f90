! This module defines the data structure 'conductor', which models the physical state of a metallic 
! conductor as a function of position and energy. The purpose of this structure is mainly to be used
! as a base class for more interesting material classes, such as those that describe superconductors
! and ferromagnets, or to be used in conjunction with such materials to describe hybrid structures.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-11
! Updated: 2015-07-16

module module_conductor
  use module_precision
  use module_spin
  use module_state
  implicit none

  ! Class declaration
  type conductor
    ! Physical parameters: these control the physical characteristics of the material (should be modified by the user before simulations)
    real(dp)                  :: thouless      = 1.000_dp ! Thouless energy of the material (ratio of the diffusion constant to the squared material length)
    real(dp)                  :: scattering    = 0.001_dp ! Imaginary energy term (this models inelastic scattering processes and stabilizes the BVP solver)

    class(conductor), pointer :: material_a               ! Material connected to this one at the left  interface (default: null pointer)
    class(conductor), pointer :: material_b               ! Material connected to this one at the right interface (default: null pointer)
    real(dp)                  :: conductance_a = 0.0_dp   ! Tunneling conductance of the left interface  (relative to the bulk conductance of this material)
    real(dp)                  :: conductance_b = 0.0_dp   ! Tunneling conductance of the right interface (relative to the bulk conductance of this material)
  
    ! Simulation parameters: these control the behaviour of the boundary value problem solver (can be modified by the user if necessary)
    integer                   :: scaling       = 128      ! How much larger than the initial position mesh the largest internal mesh can be (range: >1)
    integer                   :: information   = 0        ! Amount of debug information that the BVP solver should write to standard out (range: [-1,2])
    real(dp)                  :: tolerance     = 1e-6_dp  ! Error tolerance level (determines the maximum defect or global error, depends on the above)

    ! Core structure: these essential variables store the physical state of the material (can be accessed by the user)
    type(state), allocatable  :: state(:,:)               ! Physical state as a function of energy and position (respectively)
    real(dp),    allocatable  :: energy(:)                ! Discretized energy domain 
    real(dp),    allocatable  :: location(:)              ! Discretized position domain 


    ! Temp structure: these private variables are only used by internal subroutines (can not be accessed by the user)
    type(state), private      :: state_a                  ! Temporary storage for the left  boundary condition
    type(state), private      :: state_b                  ! Temporary storage for the right boundary condition
    complex(dp), private      :: erg                      ! Temporary storage for the current working energy

    contains
    ! Physical methods: these methods extract physical information about the material (invoked after simulations)
    !procedure :: get_dos     => conductor_get_dos
    !procedure :: get_dos     => conductor_get_dos_center

    ! Simulation methods: these methods control the simulation itself (invoked during simulations)
    procedure          :: update => conductor_update      ! Updates the state of the superconductor by calling 'solve'
    procedure, private :: solve  => conductor_solve       ! Defines the boundary value problem that describes the system

    ! Core methods: these methods are used internally
    final              :: conductor_destruct              ! Class destructor
  end type

  ! Class constructor
  interface conductor
    module procedure conductor_construct_bcs
  end interface

  ! Connects materials with an interface
  interface connect
    module procedure connect_tunneling
  end interface
contains
  pure function conductor_construct_bcs(erg, gap, eps, pts) result(this)
    ! Constructs a conductor object corresponding to a BCS superconductor with a given position and energy range
    type(conductor)                   :: this     ! Conductor object that will be constructed
    real(dp),    intent(in)           :: erg(:)   ! Discretized energy domain
    real(dp),    intent(in), optional :: eps      ! Imaginary energy contribution (models inelastic scattering)
    complex(dp), intent(in), optional :: gap      ! Superconducting order parameter
    integer,     intent(in), optional :: pts      ! Number of positions (min mesh size)

    real(dp)                          :: seps
    complex(dp)                       :: sgap
    integer                           :: spts
    integer                           :: n, m

    ! Optional argument: imaginary energy contribution from inelastic scattering
    if (present(eps)) then
      seps = eps
    else
      seps = 0.001_dp
    end if
   
    ! Optional argument: superconducting order parameter
    if (present(gap)) then
      sgap = gap
    else
      sgap = (1.0_dp,0.0_dp)
    end if

    ! Optional argument: initial mesh size (minimum number of positions per energy)
    if (present(pts)) then
      spts = pts
    else
      spts = 128
    end if

    ! Make sure pointers are initialized without a target
    nullify(this%material_a)
    nullify(this%material_b)

    ! Allocate memory (if necessary)
    if (.not. allocated(this%state)) then
      allocate(this%state(size(erg), spts))
      allocate(this%location(spts))
      allocate(this%energy(size(erg)))
    end if

    ! Fill the object fields
    this%location = [ ((real(n,kind=dp)/(spts-1)), n=0,spts-1) ]
    this%energy = erg
    forall (n = 1:size(this%energy))
      forall (m = 1:size(this%location))
        this%state(n,m) = state( cmplx(erg(n),seps,kind=dp), sgap )
      end forall
    end forall
  end function

  pure subroutine conductor_destruct(this)
    ! Define the class destructor
    type(conductor), intent(inout) :: this

    ! Deallocate memory (if necessary)
    if(allocated(this%state)) then
      deallocate(this%state)
      deallocate(this%location)
      deallocate(this%energy)
    end if
  end subroutine

  subroutine conductor_update(this)
    ! This subroutine solves the Usadel equation numerically for the given position
    ! and energy range, and updates the current estimate for the state of the system.
    use bvp_m

    class(conductor), intent(inout) :: this
    real(dp)                        :: u(32,size(this%location))
    type(bvp_sol)                   :: sol
    integer                         :: n, m

    ! Loop over quasiparticle energies
    do n=1,size(this%energy)
      ! Status information
      write (*,'(a,1x,i3,1x,a,1x,i3,1x,a,2x,f8.5)') '[',n,'/',size(this%energy),']',this%energy(n)

      ! Convert all states at this energy level to real-valued state vectors
      forall (m=1:size(this%location))
        u(:,m) = this%state(n,m)
      end forall

      ! Copy the boundary conditions and energy to temporary work variables
      this%erg = cmplx(this%energy(n)/this%thouless, this%scattering/this%thouless, kind=dp)
      if (associated(this%material_a)) then
        this%state_a = this%material_a%state(n,ubound(this%material_a%state,2))
      end if
      if (associated(this%material_b)) then
        this%state_b = this%material_b%state(n,lbound(this%material_b%state,2))
      end if

      ! Initialize the boundary value problem solver
      sol = bvp_init(32, 16, this%location, u, max_num_subintervals=size(this%location)*this%scaling)

      ! Solve the differential equation
      call this%solve(sol)

      ! Use the results to update the current state of the system
      forall (m=1:size(this%location))
        this%state(n,m) = sol%y(:,m)
      end forall
    end do

    ! Clean up
    call bvp_terminate(sol)
  end subroutine

  subroutine conductor_solve(this, sol)
    ! Solves the Usadel equation for a normal metallic conductor
    use bvp_m

    class(conductor), intent(in)    :: this
    type(bvp_sol),    intent(inout) :: sol
    
    sol = bvp_solver(sol, ode, bc, method=6, error_control=1, tol=this%tolerance, trace=this%information)
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
      d2g  = (-2.0_dp)*dg*Nt*gt*dg - (0.0_dp,2.0_dp)*this%erg*g
      d2gt = (-2.0_dp)*dgt*N*g*dgt - (0.0_dp,2.0_dp)*this%erg*gt
            
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
      if (associated(this%material_a)) then
        ! State in the material to the left
        g0   = this%state_a%g
        gt0  = this%state_a%gt
        dg0  = this%state_a%dg
        dgt0 = this%state_a%dgt

        ! Calculate the normalization matrices
        N0   = spin_inv( pauli0 - g0*gt0 )
        Nt0  = spin_inv( pauli0 - gt0*g0 )

        ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
        bca(1: 8) = dg1  - this%conductance_a*( pauli0 - g1*gt0 )*N0*(  g1  - g0  )
        bca(9:16) = dgt1 - this%conductance_a*( pauli0 - gt1*g0 )*Nt0*( gt1 - gt0 )
      else
        ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
        bca(1: 8) = dg1
        bca(9:16) = dgt1
      end if
        
      ! Right boundary condition
      if (associated(this%material_b)) then
        ! State in the material to the right
        g3   = this%state_b%g;
        gt3  = this%state_b%gt;
        dg3  = this%state_b%dg;
        dgt3 = this%state_b%dgt;

        ! Calculate the normalization matrices
        N3   = spin_inv( pauli0 - g3*gt3 );
        Nt3  = spin_inv( pauli0 - gt3*g3 );
        
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

  subroutine connect_tunneling(material_a, material_b, conductance_a, conductance_b)
    ! Creates an interface between two conductive materials.
    class(conductor), target, intent(inout) :: material_a      ! This object represents the left  material
    class(conductor), target, intent(inout) :: material_b      ! This object represents the right material
    real(dp),                 intent(in)    :: conductance_a   ! Tunneling conductance of the interface (relative to the left  bulk conductance)
    real(dp),                 intent(in)    :: conductance_b   ! Tunneling conductance of the interface (relative to the right bulk conductance)
    
    ! Update the internal material pointers
    material_a%material_b => material_b
    material_b%material_a => material_a

    ! Update the interface parameters
    material_a%conductance_b = conductance_a
    material_b%conductance_a = conductance_b
  end subroutine
end module
