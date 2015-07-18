! This module defines the data structure 'conductor', which models the physical state of a metallic 
! conductor as a function of position and energy. The purpose of this structure is mainly to be used
! as a base class for more interesting material classes, such as those that describe superconductors
! and ferromagnets, or to be used in conjunction with such materials to describe hybrid structures.
! The module also defines some functions that can be applied to any objects of the class 'conductor'.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-11
! Updated: 2015-07-17

module module_conductor
  use module_precision
  use module_spin
  use module_state
  implicit none

  ! Type declaration
  type conductor
    ! Physical parameters: these control the physical characteristics of the material (should be modified by the user)
    real(dp)                  :: thouless      = 1.000_dp                    ! Thouless energy of the material (ratio of the diffusion constant to the squared material length)
    real(dp)                  :: scattering    = 0.001_dp                    ! Imaginary energy term (this models inelastic scattering processes and stabilizes the BVP solver)
    class(conductor), pointer :: material_a                                  ! Material connected to this one at the left  interface (default: null pointer)
    class(conductor), pointer :: material_b                                  ! Material connected to this one at the right interface (default: null pointer)
    real(dp)                  :: conductance_a = 0.0_dp                      ! Tunneling conductance of the left interface  (relative to the bulk conductance of this material)
    real(dp)                  :: conductance_b = 0.0_dp                      ! Tunneling conductance of the right interface (relative to the bulk conductance of this material)
  
    ! Simulation parameters: these control the behaviour of the boundary value problem solver (can be modified by the user)
    integer                   :: scaling       = 128                         ! How much larger than the initial position mesh the largest internal mesh can be (range: >1)
    integer                   :: information   = 0                           ! Amount of status information that should be written to standard out (range: [-1,2])
    real(dp)                  :: tolerance     = 1e-6_dp                     ! Error tolerance level (determines the maximum allowed defect in the solution)

    ! Core structure: these essential variables store the physical state of the material (can be accessed by the user)
    type(state), allocatable  :: state(:,:)                                  ! Physical state as a function of energy and position
    real(dp),    allocatable  :: energy(:)                                   ! Discretized  energy  domain that will be considered
    real(dp),    allocatable  :: location(:)                                 ! Discretized position domain that will be considered

    ! Temp structure: these private variables are only used by internal subroutines (should not be accessed by the user)
    type(state)               :: state_a                                     ! Temporary storage for the left  boundary condition
    type(state)               :: state_b                                     ! Temporary storage for the right boundary condition
    complex(dp)               :: erg                                         ! Temporary storage for the current working energy

    contains
    ! Simulation methods: these methods control the simulation process (should be invoked by the user)
    procedure          :: update             => conductor_update             ! Updates the state of the material

    ! 
    procedure, private :: usadel_equation    => conductor_usadel_equation    ! Differential equation that describes the conductor
    procedure, private :: interface_vacuum_a => conductor_interface_vacuum_a ! Defines the left  boundary condition for a vacuum interface
    procedure, private :: interface_vacuum_b => conductor_interface_vacuum_b ! Defines the right boundary condition for a vacuum interface
    procedure, private :: interface_tunnel_a => conductor_interface_tunnel_a ! Defines the left  boundary condition for a tunnel interface
    procedure, private :: interface_tunnel_b => conductor_interface_tunnel_b ! Defines the right boundary condition for a tunnel interface
    procedure, private :: internals_update   => conductor_internals_update   ! 

    ! Output methods: these methods are used to export physical results to files (can be invoked by the user)
    procedure          :: write_dos          => conductor_write_dos          ! Writes the density of states to a given output unit

    ! Core methods: these methods are only used by internal subroutines (should not be directly invoked by the user)
    final              :: conductor_destruct                                 ! Destructor that deallocates dynamic memory
  end type

  ! Constructs a conductor and initializes it a superconducting state
  interface conductor
    module procedure conductor_construct_bcs
  end interface

  ! Constructs an interface by connecting two conductor objects
  interface connect
    module procedure connect_tunneling
  end interface
contains
  pure function conductor_construct_bcs(energy, gap, scattering, points) result(this)
    ! Constructs a conductor object corresponding to a BCS superconductor with a given position and energy range
    type(conductor)                   :: this        ! Conductor object that will be constructed
    real(dp),    intent(in)           :: energy(:)   ! Discretized energy domain that will be used
    real(dp),    intent(in), optional :: scattering  ! Imaginary energy term (default: 0.001)
    complex(dp), intent(in), optional :: gap         ! Superconducting gap (default: 1.0)
    integer,     intent(in), optional :: points      ! Number of positions (default: 128)

    real(dp)                          :: seps
    complex(dp)                       :: sgap
    integer                           :: spts
    integer                           :: n, m

    ! Optional argument: imaginary energy contribution due to inelastic scattering
    if (present(scattering)) then
      seps = scattering
    else
      seps = 0.001_dp
    end if
   
    ! Optional argument: superconducting order parameter
    if (present(gap)) then
      sgap = gap
    else
      sgap = (1.0_dp,0.0_dp)
    end if

    ! Optional argument: number of positions to use
    if (present(points)) then
      spts = points
    else
      spts = 128
    end if

    ! Make sure pointers are initialized to null pointers
    nullify(this%material_a)
    nullify(this%material_b)

    ! Allocate memory (if necessary)
    if (.not. allocated(this%state)) then
      allocate(this%state(size(energy), spts))
      allocate(this%location(spts))
      allocate(this%energy(size(energy)))
    end if

    ! Fill the object fields
    this%location = [ ((real(n,kind=dp)/real(spts-1,kind=dp)), n=0,spts-1) ]
    this%energy = energy
    forall (n = 1:size(this%energy))
      forall (m = 1:size(this%location))
        this%state(n,m) = state( cmplx(energy(n),seps,kind=dp), sgap )
      end forall
    end forall
  end function

  pure subroutine conductor_destruct(this)
    ! Define the type destructor
    type(conductor), intent(inout) :: this

    ! Deallocate memory (if necessary)
    if(allocated(this%state)) then
      deallocate(this%state)
      deallocate(this%location)
      deallocate(this%energy)
    end if
  end subroutine

  subroutine conductor_update(this)
    ! This subroutine updates the current estimate for the state of the system.
    use bvp_m

    class(conductor), intent(inout) :: this                       ! Object that will be updated
    real(dp)                        :: u(32,size(this%location))  ! Real state vector required by the BVP solver
    type(bvp_sol)                   :: sol                        ! Object with information about the BVP solution
    integer                         :: n, m                       ! Internal loop variables

    do n=1,size(this%energy)
      ! Status information
      if (this%information >= 0) then
        write (*,'(3x,a,1x,i3,1x,a,1x,i3,1x,a,1x,f8.5)') '[',n,'/',size(this%energy),']',this%energy(n)
      end if

      ! Convert all states at this energy level to real-valued state vectors
      forall (m=1:size(this%location))
        u(:,m) = this%state(n,m)
      end forall

      ! Copy the complex energy and boundary conditions to internal work variables
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
      sol = bvp_solver(sol, ode, bc, method=6, error_control=1, tol=this%tolerance, trace=this%information)

      ! Use the results to update the current state of the system
      forall (m=1:size(this%location))
        this%state(n,m) = sol%y(:,m)
      end forall

      ! Update other internal variables if necessary
      call this%internals_update
    end do

    ! Clean up
    call bvp_terminate(sol)
  contains
    subroutine ode(z, u, f)
      ! Definition of the differential equation u'=f(z,u)
      real(dp), intent(in)  :: z
      real(dp), intent(in)  :: u(32)
      real(dp), intent(out) :: f(32)

      type(spin)            :: g, gt, dg, dgt, d2g, d2gt

      ! Extract the Riccati parameters
      g   = u( 1: 8)
      gt  = u( 9:16)
      dg  = u(17:24)
      dgt = u(25:32)

      ! Calculate the second derivatives of the Riccati parameters
      call this%usadel_equation(z, g, gt, dg, dgt, d2g, d2gt)
            
      ! Pack the results into a state vector
      f( 1: 8) = dg
      f( 9:16) = dgt
      f(17:24) = d2g
      f(25:32) = d2gt
    end subroutine

    subroutine bc(ua, ub, bca, bcb)
      ! Definition of the boundary conditions bca=g(ua) and bcb=g(ub)
      real(dp), intent(in)  :: ua(32)
      real(dp), intent(in)  :: ub(32)
      real(dp), intent(out) :: bca(16)
      real(dp), intent(out) :: bcb(16)

      type(spin)            :: g1, gt1, dg1, dgt1, r1, rt1
      type(spin)            :: g2, gt2, dg2, dgt2, r2, rt2

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

      ! Calculate residuals from the boundary conditions
      if (associated(this%material_a)) then
        call this%interface_tunnel_a(g1, gt1, dg1, dgt1, r1, rt1)
      else
        call this%interface_vacuum_a(g1, gt1, dg1, dgt1, r1, rt1)
      end if

      if (associated(this%material_b)) then
        call this%interface_tunnel_b(g2, gt2, dg2, dgt2, r2, rt2)
      else
        call this%interface_vacuum_b(g2, gt2, dg2, dgt2, r2, rt2)
      end if

      ! Pack the results into state vectors
      bca(1: 8) = r1
      bca(9:16) = rt1
      bcb(1: 8) = r2
      bcb(9:16) = rt2
    end subroutine
  end subroutine

  subroutine conductor_usadel_equation(this, z, g, gt, dg, dgt, d2g, d2gt)
    ! Use the Usadel equation to calculate the second derivatives of the Riccati parameters at point z.
    class(conductor), intent(in)  :: this
    real(dp),         intent(in)  :: z
    type(spin),       intent(in)  :: g, gt, dg, dgt
    type(spin),       intent(out) :: d2g, d2gt
    type(spin)                    :: N, Nt

    ! Calculate the normalization matrices
    N   = spin_inv( pauli0 - g*gt )
    Nt  = spin_inv( pauli0 - gt*g )

    ! Calculate the second derivatives of the Riccati parameters
    d2g  = (-2.0_dp)*dg*Nt*gt*dg - (0.0_dp,2.0_dp)*this%erg*g
    d2gt = (-2.0_dp)*dgt*N*g*dgt - (0.0_dp,2.0_dp)*this%erg*gt
  end subroutine

  subroutine conductor_interface_vacuum_a(this, g1, gt1, dg1, dgt1, r1, rt1)
    ! Defines a vacuum boundary condition for the left interface
    class(conductor), intent(in)  :: this
    type(spin),       intent(in)  :: g1, gt1, dg1, dgt1
    type(spin),       intent(out) :: r1, rt1

    r1  = dg1
    rt1 = dgt1
  end subroutine

  subroutine conductor_interface_vacuum_b(this, g2, gt2, dg2, dgt2, r2, rt2)
    ! Defines a vacuum boundary condition for the right interface
    class(conductor), intent(in)  :: this
    type(spin),       intent(in)  :: g2, gt2, dg2, dgt2
    type(spin),       intent(out) :: r2, rt2

    r2  = dg2
    rt2 = dgt2
  end subroutine

  subroutine conductor_interface_tunnel_a(this, g1, gt1, dg1, dgt1, r1, rt1)
    ! Defines a tunneling boundary condition for the left interface 
    class(conductor), intent(in)  :: this
    type(spin),       intent(out) :: r1, rt1
    type(spin),       intent(in)  :: g1, gt1, dg1, dgt1
    type(spin)                    :: g0, gt0, dg0, dgt0, N0, Nt0

    ! Rename the state in the material to the left
    g0   = this%state_a%g
    gt0  = this%state_a%gt
    dg0  = this%state_a%dg
    dgt0 = this%state_a%dgt

    ! Calculate the normalization matrices
    N0   = spin_inv( pauli0 - g0*gt0 )
    Nt0  = spin_inv( pauli0 - gt0*g0 )

    ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
    r1  = dg1  - this%conductance_a*( pauli0 - g1*gt0 )*N0*(  g1  - g0  )
    rt1 = dgt1 - this%conductance_a*( pauli0 - gt1*g0 )*Nt0*( gt1 - gt0 )
  end subroutine

  subroutine conductor_interface_tunnel_b(this, g2, gt2, dg2, dgt2, r2, rt2)
    ! Defines a tunneling boundary condition for the right interface 
    class(conductor), intent(in)  :: this
    type(spin),       intent(out) :: r2, rt2
    type(spin),       intent(in)  :: g2, gt2, dg2, dgt2
    type(spin)                    :: g3, gt3, dg3, dgt3, N3, Nt3

    ! Rename the state in the material to the right
    g3   = this%state_b%g;
    gt3  = this%state_b%gt;
    dg3  = this%state_b%dg;
    dgt3 = this%state_b%dgt;
  
    ! Calculate the normalization matrices
    N3   = spin_inv( pauli0 - g3*gt3 );
    Nt3  = spin_inv( pauli0 - gt3*g3 );

    ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
    r2  = dg2  - this%conductance_b*( pauli0 - g2*gt3 )*N3*(  g3  - g2  )
    rt2 = dgt2 - this%conductance_b*( pauli0 - gt2*g3 )*Nt3*( gt3 - gt2 )
  end subroutine

  subroutine conductor_internals_update(this)
    class(conductor), intent(inout) :: this

    continue
  end subroutine

  subroutine connect_tunneling(material_a, material_b, conductance_a, conductance_b)
    ! Creates a tunneling interface between two conductive materials.
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

  subroutine conductor_write_dos(this, unit)
    ! Writes the density of states as a function of position and energy to a given output unit.
    class(conductor),   intent(in) :: this      ! Material that the density of states will be calculated from
    integer,            intent(in) :: unit      ! Output unit that determines where the information will be written
    integer                        :: n, m      ! Temporary loop variables

    ! Write the information to output
    do m=1,size(this%location)
      if (minval(this%energy) < 0.0_dp) then
        ! If we have data for both positive and negative energies, simply write out the data
        write(unit,*) (this%state(n,m)%get_dos(), n=1,size(this%energy))
      else
        ! If we only have data for positive energies, assume that the regions are symmetric
        write(unit,*) (this%state(n,m)%get_dos(), n=size(this%energy),1,-1),&
                      (this%state(n,m)%get_dos(), n=1,size(this%energy),+1)
      end if
    end do
  end subroutine
end module
