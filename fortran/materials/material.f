!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This module defines the data type 'material',  which models the state of a physical material for a discretized range of
!> positions and energies. This is an abstract type, meaning that it is not intended to be instantiated on its own, but is
!> intended as a base type for physical materials like conductors, superconductors, and ferromagnets. In other words, this
!> type defines the essential data structures and program structure, while the derived subtypes will define actual physics.
!>
!> @TODO
!>   It might be cleaner to refactor this class by moving interface_vacuum and interface_transparent here, since these
!>   are more a general description of materials than something specific to conductors. The Kupryanov-Lukichev boundary
!>   conditions should however stay in the conductor class, since these describe "physics and mot math". Furthermore, 
!>   it may be useful to move the diffusion_equation terms of the kind -2·dg·Nt·gt·dg to this class, since they're 
!>   technically a part of the Riccati parametrization itself and not something describing any particular material.
!>   With these modifications, the class can be made concrete instead of abstract, removing the need for interfaces.

module material_m
  use :: stdio_m
  use :: math_m
  use :: spin_m
  use :: propagator_m
  private

  ! Public declarations
  public :: material_conf, material_load

  ! Type declarations
  type, public, abstract :: material
    ! Material properties that affect the equilibrium state
    real(wp)                                  :: length                =  1.00_wp              !! Material length (L/ξ)
    real(wp)                                  :: thouless              =  1.00_wp              !! Thouless energy (ħD/L²)
    real(wp)                                  :: scattering            =  0.01_wp              !! Inelastic scattering (η/Δ₀)
    logical                                   :: transparent_a         =  .false.              !! Interface transparency (left)
    logical                                   :: transparent_b         =  .false.              !! Interface transparency (right)

    ! Material properties that affect the nonequilibrium state
    logical                                   :: equilibrium           =  .true.               !! Equilibrium?
    logical                                   :: transverse            =  .false.              !! Transverse potential gradients?
    real(wp)                                  :: voltage               =  0.00_wp              !! Voltage (eV/Δ₀)
    real(wp)                                  :: temperature           =  0.01_wp              !! Temperature (T/Tc)
    real(wp)                                  :: spinvoltage           =  0.00_wp              !! Spin-voltage (eVs/Δ₀)
    real(wp)                                  :: spintemperature       =  0.00_wp              !! Spin-temperature (Ts/Tc)

    ! Physical state modelled using quasiclassical propagators
    real(wp),                     allocatable :: energy(:)                                     !! Energy domain
    real(wp),                     allocatable :: location(:)                                   !! Position domain
    type(propagator),             allocatable :: propagator(:,:)                               !! Propagator values
    type(propagator),             allocatable :: backup(:,:)                                   !! Propagator backups

    ! Physical observables derived from the propagators above
    real(wp),                     allocatable :: density(:,:,:)                                !! Spin-resolved density of states
    real(wp),                     allocatable :: supercurrent(:,:)                             !! Charge, spin, heat, and spin-heat supercurrents
    real(wp),                     allocatable :: lossycurrent(:,:)                             !! Charge, spin, heat, and spin-heat dissipative currents
    real(wp),                     allocatable :: accumulation(:,:)                             !! Charge, spin, heat, and spin-heat accumulation
    complex(wp),                  allocatable :: correlation(:)                                !! Superconducting pair-correlations

    ! Hybrid structures are modeled as double-linked lists
    integer                                   :: order           =  1                          !! Simulation priority of this material
    class(material),                  pointer :: material_a      => null()                     !! Material to the left  (default: vacuum)
    class(material),                  pointer :: material_b      => null()                     !! Material to the right (default: vacuum)

    ! The package bvp_solver is used to handle differential equations, and will be controlled by the following parameters
    integer                                   :: selfconsistency =  2                          !! Selfconsistency scheme (0 = none, 1 = fixpoint iteration, 2 = Steffensen's method)
    integer                                   :: scaling         =  128                        !! Maximal mesh increase (range: 2^N, N>1)
    integer                                   :: method          =  4                          !! Runge—Kutta order (range: 2, 4, 6)
    integer                                   :: control         =  2                          !! Error control (1: defect, 2: global error, 3: 1 then 2, 4: 1 and 2)
    real(wp)                                  :: tolerance       =  1e-10_wp                   !! Error tolerance (maximum defect or global error)
    integer                                   :: information     =  0                          !! Debug information (range: [-1,2])
    real(wp)                                  :: difference      =  1e+10_wp                   !! Difference between iterations

    ! The following variables are used for input/output purposes, and should be modified by class(material) constructors
    character(len=128)                        :: type_string     =  'MATERIAL'                 !! Name of this material
  contains
    ! These methods define how to update the physical state of the material
    procedure(init),                 deferred :: init                                          !! Initializes  propagators
    procedure                                 :: update           => material_update           !! Calculates propagators
    procedure                                 :: update_diffusion => material_update_diffusion !! Calculates propagators (in equilibrium)
    procedure                                 :: update_kinetic   => material_update_kinetic   !! Calculates propagators (nonequilibrium)
    procedure(update),               deferred :: update_prehook                                !! Executed before update
    procedure(update),               deferred :: update_posthook                               !! Executed after  update

    ! These methods define the physical equations used by the update methods
    procedure(diffusion_equation),   deferred :: diffusion_equation                            !! Diffusion equation
    procedure(diffusion_equation_a), deferred :: diffusion_equation_a                          !! Boundary condition (left)
    procedure(diffusion_equation_b), deferred :: diffusion_equation_b                          !! Boundary condition (right)

    ! These methods define miscellaneous utility functions
    procedure                                 :: conf            => material_conf              !! Configures material parameters
    procedure                                 :: save            => material_save              !! Saves the state of the material
    procedure                                 :: load            => material_load              !! Loads the state of the material
  end type

  ! Interface declarations
  abstract interface
    pure subroutine init(this, gap)
      !! This interface is used for the deferred procedure init.
      import material, wp

      class(material),       intent(inout) :: this
      complex(wp), optional, intent(in)    :: gap
    end subroutine
  end interface

  abstract interface
    impure subroutine update(this)
      !! This interface is used for the deferred procedures update_prehook and update_posthook.
      import material

      class(material), intent(inout) :: this
    end subroutine
  end interface

  abstract interface
    pure subroutine diffusion_equation(this, e, z, g, gt, dg, dgt, d2g, d2gt)
      !! This interface is used for the deferred procedure diffusion_equation.
      import material, spin, wp

      class(material), intent(in)    :: this
      type(spin),      intent(in)    :: g, gt, dg, dgt
      type(spin),      intent(inout) :: d2g, d2gt
      complex(wp),     intent(in)    :: e
      real(wp),        intent(in)    :: z
    end subroutine
  end interface

  abstract interface
    pure subroutine diffusion_equation_a(this, p, a, r, rt)
      !! This interface is used for the deferred procedure diffusion_equation_a.
      import material, propagator, spin, wp

      class(material),          intent(in)    :: this
      type(propagator),         intent(in)    :: p, a
      type(spin),               intent(inout) :: r, rt
    end subroutine
  end interface

  abstract interface
    pure subroutine diffusion_equation_b(this, p, b, r, rt)
      !! This interface is used for the deferred procedure diffusion_equation_b.
      import material, propagator, spin, wp

      class(material),          intent(in)    :: this
      type(propagator),         intent(in)    :: p, b
      type(spin),               intent(inout) :: r, rt
    end subroutine
  end interface
contains

  !--------------------------------------------------------------------------------!
  !                    IMPLEMENTATION OF STATE UPDATE METHODS                      !
  !--------------------------------------------------------------------------------!

  impure subroutine material_update(this, bootstrap)
    !! This subroutine updates the state of the material by solving diffusion equations for the equilibrium
    !! propagators, kinetic equations for the nonequilibrium propagators, and then calculating observables.
    class(material),   intent(inout) :: this                       !! Material that will be updated
    logical, optional, intent(in)    :: bootstrap                  !! Disables calculation of observables

    ! Only update materials with a positive order
    if (this % order <= 0) then
      return
    end if

    ! Call the prehook method
    call this % update_prehook

    ! Status information
    if (this % information >= 0) then
      block
        character(len=10) :: str
        write(str, '(2x,"[1m#",i0)') this % order
        write(stdout,'(a,a,20x)') str, trim(this % type_string)
      end block
    end if

    ! Reset the difference since last update to zero
    this % difference = 0.0_wp

    ! Solve the diffusion equations
    call this % update_diffusion

    ! Solve the kinetic equations
    if (.not. this % equilibrium) then
      call this % update_kinetic
    end if

    ! Status information
    if (this % information >= 0) then
      write(stdout,'(6x,a,f10.8,a)') 'Max change: ',this % difference,'                                        '
      flush(stdout)
    end if

    ! Stop here if bootstrapping
    if (present(bootstrap)) then
      if (bootstrap) then
        return
      end if
    end if

    ! Call the posthook method
    call this % update_posthook
  end subroutine

  impure subroutine material_update_diffusion(this)
    !! This subroutine calculates the equilibrium propagators of the material from the diffusion equation.
    use bvp_m

    class(material),              intent(inout) :: this
    type(bvp_sol)                               :: sol
    real(wp), dimension(32,size(this%location)) :: u, d
    type(propagator)                            :: a, b
    integer                                     :: n, m
    complex(wp)                                 :: e

    ! Loop over energies
    do n=size(this%energy),1,-1
      ! Status information
      if (this%information >= 0) then
        write(stdout,'(6x,a,1x,i4,1x,a,1x,i4,1x,a,f0.5,a1)',advance='no') &
          '[',n,'/',size(this%energy),']  E = ',this%energy(n), achar(13)
        flush(stdout)
      end if

      ! Construct the state vector from the Riccati parameters and their derivatives
      do m=1,size(this%location)
        call pack(this % propagator(n,m), d(:,m))
      end do

      ! If the difference between iterations is small, use the results from the previous
      ! iteration as an initial guess. If not, use the results form the previous energy.
      if (this % difference < 0.05_wp) then
        u = d
      end if

      ! Calculate the complex normalized energy
      e = cx(this%energy(n), this%scattering)/this%thouless

      ! Update the boundary conditions (left)
      a = propagator()
      if (associated(this % material_a)) then
        associate(other => this % material_a % propagator)
          a = other(n,ubound(other,2))
        end associate
      end if

      ! Update the boundary conditions (right)
      b = propagator()
      if (associated(this % material_b)) then
        associate(other => this % material_b % propagator)
          b = other(n,lbound(other,2))
        end associate
      end if

      ! Initialize bvp_solver
      sol = bvp_init(32, 16, this%location, u, max_num_subintervals=(size(this%location)*this%scaling))

      ! Solve the differential equation
      sol = bvp_solver(sol, ode, bc, method=this%method, error_control=this%control, &
                       tol=this%tolerance, trace=this%information, stop_on_fail=.false.)

      ! Check if the calculation succeeded
      if (sol % info /= 0) then
        call error('Failed to converge! This is usually because of an ill-posed problem.')
      end if

      ! Use the results to update the state
      call bvp_eval(sol, this % location, u)
      do m=1,size(this%location)
        call unpack(this % propagator(n,m), u(:,m))
      end do

      ! Update the difference vector
      d = abs(u - d)

      ! Update the maximal difference since last update
      this % difference = max(this % difference, maxval(d))

      ! Clean up after bvp_solver
      call bvp_terminate(sol)
    end do
  contains
    pure subroutine ode(z, u, f)
      !! Definition of the differential equation du/dz=f(z,u).
      real(wp), intent(in)  :: z
      real(wp), intent(in)  :: u(32)
      real(wp), intent(out) :: f(32)
      type(spin)            :: g, gt, dg, dgt, d2g, d2gt

      ! Extract the Riccati parameters
      g   = u( 1: 8)
      gt  = u( 9:16)
      dg  = u(17:24)
      dgt = u(25:32)

      ! Calculate the second-derivatives of the Riccati parameters
      call this%diffusion_equation(e, z, g, gt, dg, dgt, d2g, d2gt)

      ! Pack the results into a state vector
      f( 1: 8) = dg
      f( 9:16) = dgt
      f(17:24) = d2g
      f(25:32) = d2gt
    end subroutine

    pure subroutine bc(ua, ub, bca, bcb)
      !! Definition of the boundary conditions bca=g(ua) and bcb=g(ub).
      real(wp), intent(in)  :: ua(32)
      real(wp), intent(in)  :: ub(32)
      real(wp), intent(out) :: bca(16)
      real(wp), intent(out) :: bcb(16)

      type(propagator)      :: pa, pb
      type(spin)            :: ra, rta
      type(spin)            :: rb, rtb

      ! State at the left end of the material
      call unpack(pa, ua)

      ! State at the right end of the material
      call unpack(pb, ub)

      ! Calculate residuals from the boundary conditions (left)
      if (this % transparent_a) then
        ! Transparent interface
        ra  = pa % g  - a % g
        rta = pa % gt - a % gt
      else
        ! Customized interface
        call this % diffusion_equation_a(pa, a, ra, rta)
      end if

      ! Calculate residuals from the boundary conditions (right)
      if (this % transparent_b) then
        ! Transparent interface
        rb  = pb % g  - b % g
        rtb = pb % gt - b % gt
      else
        ! Customized interface
        call this % diffusion_equation_b(pb, b, rb, rtb)
      end if

      ! Pack the results into state vectors
      bca(1: 8) = ra
      bca(9:16) = rta
      bcb(1: 8) = rb
      bcb(9:16) = rtb
    end subroutine

    pure subroutine pack(a, b)
      !! Defines assignment from a propagator object to a real vector.
      type(propagator),        intent(in)  :: a
      real(wp), dimension(32), intent(out) :: b

      b( 1: 8) = a % g
      b( 9:16) = a % gt
      b(17:24) = a % dg
      b(25:32) = a % dgt
    end subroutine

    pure subroutine unpack(a, b)
      !! Defines assignment from a real vector to a propagator object.
      type(propagator),        intent(inout) :: a
      real(wp), dimension(32), intent(in)    :: b

      ! Unpack the vector elements
      a % g   = b( 1: 8) 
      a % gt  = b( 9:16) 
      a % dg  = b(17:24) 
      a % dgt = b(25:32) 

      ! Update normalization matrices
      a % N  = inverse( pauli0 - a%g  * a%gt )
      a % Nt = inverse( pauli0 - a%gt * a%g  )
    end subroutine
  end subroutine

  impure subroutine material_update_kinetic(this)
    !! This subroutine calculates the nonequilibrium propagators of the material from the kinetic equation.
    use bvp_m

    ! Material object
    class(material),                 intent(inout) :: this

    ! Numerical solver
    type(bvp_sol)                                  :: sol

    ! State vectors
    real(wp), dimension(16,size(this%location))    :: u, d
    real(wp), dimension(16)                        :: a, b

    ! Jacobian matrices
    real(wp), dimension(16,16,size(this%location)) :: ju
    real(wp), dimension(16,16)                     :: ja, jb

    ! Loop variables
    integer                                        :: n, m

    ! Loop over energies
    do n=size(this%energy),1,-1
      ! Status information
      if (this%information >= 0) then
        write(stdout,'(6x,a,1x,i4,1x,a,1x,i4,1x,a,f0.5,a1)',advance='no') &
          '[',n,'/',size(this%energy),']  E = ',this%energy(n), achar(13)
        flush(stdout)
      end if

      ! Construct the state vector from the distribution and its derivative
      do m=1,size(this%location)
        call pack(this % propagator(n,m), d(:,m))
      end do

      ! Use the results from the previous iteration as an initial guess
      u = d

      ! Update the Jacobian matrix at each position
      do m=1,size(this%location)
        ! Trivial:  ∂(H)'/∂H
        ju(:8,:8,m) = 0

        ! Trivial:  ∂(H)'/∂H'
        ju(:8,9:,m) = identity(8)

        ! Physical: ∂(H')'/∂H
        ju(9:,:8,m) = 0

        ! Physical: ∂(H')'/∂H'
        ju(9:,9:,m) = 0
      end do

      ! Update the boundary conditions (left)
      a = 0
      if (associated(this % material_a)) then
        associate( other => this % material_a % propagator )
          call pack(other(n,ubound(other,2)), a)
        end associate
      end if

      ! Update the boundary conditions (right)
      b = 0
      if (associated(this % material_b)) then
        associate( other => this % material_b % propagator )
          call pack(other(n,lbound(other,2)), b)
        end associate
      end if

      ! Initialize bvp_solver
      sol = bvp_init(16, 8, this%location, u, max_num_subintervals=(size(this%location)*this%scaling))

      ! Solve the differential equation
      sol = bvp_solver(sol, ode, bc, dfdy=jac, method=this%method, error_control=this%control, &
                       tol=this%tolerance, trace=this%information, stop_on_fail=.false.)

      ! Check if the calculation succeeded
      if (sol % info /= 0) then
        call error('Failed to converge! This is usually because of an ill-posed problem.')
      end if

      ! Use the results to update the state
      call bvp_eval(sol, this % location, u)
      do m=1,size(this%location)
        call unpack(this % propagator(n,m), u(:,m))
      end do

      ! Update the difference vector
      d = abs(u - d)

      ! Update the maximal difference since last update
      this % difference = max(this % difference, maxval(d))

      ! Clean up after bvp_solver
      call bvp_terminate(sol)
    end do
  contains
    pure subroutine ode(z, u, f)
      !! Definition of the differential equation du/dz=f(z,u).
      real(wp),                intent(in)  :: z
      real(wp), dimension(16), intent(in)  :: u
      real(wp), dimension(16), intent(out) :: f
      real(wp), dimension(16,16)           :: j

      ! Calculate the Jacobian matrix
      call jac(z, u, j)

      ! Calculate the state derivative
      f = matmul(j, u)
    end subroutine

    pure subroutine jac(z, u, j)
      !! Jacobian of the differential equation j=∂f/∂u.
      !! @TODO: Calculate M, Q, K to find the true Jacobian.
      !! @TODO: Call the kinetic equation to find self-energies.
      real(wp),                   intent(in)  :: z
      real(wp), dimension(16),    intent(in)  :: u
      real(wp), dimension(16,16), intent(out) :: j
      integer                                 :: m, n

      ! Calculate the index corresponding to the given location
      m = size(this % location) - 1
      n = floor(z*m + 1)

      ! Interpolate the Jacobian at this position
      if (n <= 1) then
        j(:,:) = ju(:,:,1)
      else
        j(:,:) = ju(:,:,n-1) + (ju(:,:,n)-ju(:,:,n-1))*(z*m-(n-2))
      end if
    end subroutine

    pure subroutine bc(ua, ub, bca, bcb)
      !! Definition of the boundary conditions bca=g(ua) and bcb=g(ub).
      !! @TODO: Implement support for general boundary conditions.
      !! @TODO: Call kinetic_equation_a and kinetic_equation_b.
      !! @TODO: Calculate the Jacobian of the boundary conditions?
      real(wp), intent(in)  :: ua(16)
      real(wp), intent(in)  :: ub(16)
      real(wp), intent(out) :: bca(8)
      real(wp), intent(out) :: bcb(8)

      ! Calculate residuals from the boundary conditions (left)
      if (this % transparent_a) then
        ! Transparent interface
        bca(1:8) = ua(1:8) - a(1:8)
      ! else
      ! ! Customized interface
      !   call this % kinetic_equation_a(ua, a)
      end if

      if (this % transparent_b) then
        ! Transparent interface
        bcb(1:8) = ub(1:8) - b(1:8)
      ! else
      ! ! Customized interface
      !   call this % kinetic_equation_b(ub, b)
      end if
    end subroutine

    pure subroutine pack(a, b)
      !! Defines assignment from a propagator object to a real vector.
      type(propagator),        intent(in)  :: a
      real(wp), dimension(16), intent(out) :: b

      b( 1: 8) = a % h
      b( 9:16) = a % dh
    end subroutine

    pure subroutine unpack(a, b)
      !! Defines assignment from a real vector to a propagator object.
      type(propagator),        intent(inout) :: a
      real(wp), dimension(16), intent(in)    :: b

      a % h   = b( 1: 8) 
      a % dh  = b( 9:16) 
    end subroutine
  end subroutine

  !--------------------------------------------------------------------------------!
  !                      IMPLEMENTATION OF UTILITY METHODS                         !
  !--------------------------------------------------------------------------------!

  impure subroutine material_conf(this, key, val)
    !! Configure a material property based on a key-value pair.
    use :: evaluate_m

    class(material), intent(inout) :: this
    character(*),    intent(in)    :: key
    character(*),    intent(in)    :: val

    select case(key)
      case("length")
        call evaluate(val, this%length)
        this%thouless = 1/(this%length**2 + eps)
      case("scattering")
        call evaluate(val, this%scattering)
      case("temperature")
        call evaluate(val, this%temperature)
        call this % init()
      case("voltage")
        call evaluate(val, this%voltage)
        call this % init()
      case("spinvoltage")
        call evaluate(val, this%spinvoltage)
        call this % init()
      case("spintemperature")
        call evaluate(val, this%spintemperature)
        call this % init()
      case("transverse")
        call evaluate(val, this%transverse)
        call this % init()
      case("transparent_a")
        call evaluate(val, this%transparent_a)
      case("transparent_b")
        call evaluate(val, this%transparent_b)
      case("order")
        call evaluate(val, this%order)
        if (this%order > 16) then
          call error("The order of the material must be be maximum 16.")
        else if (this%order < 0) then
          call error("The order of the material must be be minimum 0.")
        end if
      case("selfconsistency")
        call evaluate(val, this%selfconsistency)
        if (this % selfconsistency < 0 .or. this % selfconsistency > 2) then
          call error("The selfconsistency scheme should be in the range [0,2].")
        end if
      case("equilibrium")
        call evaluate(val, this%equilibrium)
      case default
        call warning("Unknown material option '" // key // "' ignored.")
    end select
  end subroutine

  impure subroutine material_save(this)
    !! Save a backup of the current material state.
    class(material), intent(inout) :: this

    ! Make sure the backup exists
    if (.not. allocated(this%backup)) then
      allocate(this%backup(size(this%propagator,1),size(this%propagator,2)))
    end if

    ! Make a backup of the propagators
    call this % propagator % save(this % backup)
  end subroutine

  impure subroutine material_load(this)
    !! Load a backup of a previous material state.
    class(material), intent(inout) :: this
    integer                        :: info

    ! Load the saved propagators
    if (allocated(this%backup)) then
      call this % propagator % load(this % backup)
    end if

    ! Disable status messages
    info = this%information
    if (this%information >= 0) then
      this%information = -1
    end if

    ! Silently recalculate derived properties
    call this%update_posthook

    ! Reenable status messages
    this%information = info
  end subroutine
end module
