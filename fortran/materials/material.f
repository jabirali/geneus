!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This module defines the data type 'material',  which models the state of a physical material for a discretized range of
!> positions and energies. This is an abstract type, meaning that it is not intended to be instantiated on its own, but is
!> intended as a base type for physical materials like conductors, superconductors, and ferromagnets. In other words, this
!> type defines the essential data structures and program structure, while the derived subtypes will define actual physics.

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
    ! Physical properties of a diffusive material
    real(wp)                                  :: length                =  1.00_wp           !! Material length (L/ξ)
    real(wp)                                  :: thouless              =  1.00_wp           !! Thouless energy (ħD/L²)
    real(wp)                                  :: temperature           =  0.00_wp           !! Material temperature (T/Tc)
    real(wp)                                  :: scattering            =  0.01_wp           !! Inelastic scattering (η/Δ₀)
    real(wp)                                  :: conductance_a         =  0.00_wp           !! Interface conductance (left)
    real(wp)                                  :: conductance_b         =  0.00_wp           !! Interface conductance (right)

    ! Physical state modelled using quasiclassical propagators
    real(wp),                     allocatable :: energy(:)                                  !! Energy domain
    real(wp),                     allocatable :: location(:)                                !! Position domain
    type(propagator),             allocatable :: propagator(:,:)                            !! Propagator values
    type(propagator),             allocatable :: backup(:,:)                                !! Propagator backups

    ! Physical observables derived from the propagators above
    real(wp),                     allocatable :: density(:,:,:)                             !! Spin-resolved density of states
    real(wp),                     allocatable :: supercurrent(:,:)                          !! Charge, spin, heat, and spin-heat supercurrents
    real(wp),                     allocatable :: lossycurrent(:,:)                          !! Charge, spin, heat, and spin-heat dissipative currents
    real(wp),                     allocatable :: accumulation(:,:)                          !! Charge, spin, heat, and spin-heat accumulation
    complex(wp),                  allocatable :: correlation(:)                             !! Superconducting pair-correlations

    ! Hybrid structures are modeled as double-linked lists
    integer                                   :: order           =  1                       !! Simulation priority of this material
    class(material),                  pointer :: material_a      => null()                  !! Material to the left  (default: vacuum)
    class(material),                  pointer :: material_b      => null()                  !! Material to the right (default: vacuum)

    ! The package bvp_solver is used to handle differential equations, and will be controlled by the following parameters
    integer                                   :: scaling         =  128                     !! Maximal mesh increase (range: 2^N, N>1)
    integer                                   :: method          =  4                       !! Runge—Kutta order (range: 2, 4, 6)
    integer                                   :: control         =  2                       !! Error control (1: defect, 2: global error, 3: 1 then 2, 4: 1 and 2)
    real(wp)                                  :: tolerance       =  1e-10_wp                !! Error tolerance (maximum defect or global error)
    integer                                   :: information     =  0                       !! Debug information (range: [-1,2])
    real(wp)                                  :: difference      =  1e+10_wp                !! Difference between iterations

    ! The following variables are used for input/output purposes, and should be modified by class(material) constructors
    character(len=128)                        :: type_string     =  'MATERIAL'              !! Name of this material
  contains
    ! These methods define how to update the physical state of the material
    procedure(init),                 deferred :: init                                       !! Initializes  propagators
    procedure                                 :: update          => material_update         !! Recalculates propagators
    procedure(update),               deferred :: update_prehook                             !! Executed before update
    procedure(update),               deferred :: update_posthook                            !! Executed after  update

    ! These methods define the physical equations used by the update methods
    procedure(diffusion_equation),   deferred :: diffusion_equation                         !! Diffusion equation
    procedure(interface_equation_a), deferred :: interface_equation_a                       !! Boundary condition (left)
    procedure(interface_equation_b), deferred :: interface_equation_b                       !! Boundary condition (right)

    ! These methods define miscellaneous utility functions
    procedure                                 :: conf            => material_conf           !! Configures material parameters
    procedure                                 :: save            => material_save           !! Saves the state of the material
    procedure                                 :: load            => material_load           !! Loads the state of the material
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
    pure subroutine interface_equation_a(this, a, g, gt, dg, dgt, r, rt)
      !! This interface is used for the deferred procedure interface_equation_a.
      import material, propagator, spin, wp

      class(material),          intent(in)    :: this
      type(propagator),         intent(in)    :: a
      type(spin),               intent(in)    :: g, gt, dg, dgt
      type(spin),               intent(inout) :: r, rt
    end subroutine
  end interface

  abstract interface
    pure subroutine interface_equation_b(this, b, g, gt, dg, dgt, r, rt)
      !! This interface is used for the deferred procedure interface_equation_b.
      import material, propagator, spin, wp

      class(material),          intent(in)    :: this
      type(propagator),         intent(in)    :: b
      type(spin),               intent(in)    :: g, gt, dg, dgt
      type(spin),               intent(inout) :: r, rt
    end subroutine
  end interface
contains

  !--------------------------------------------------------------------------------!
  !                    IMPLEMENTATION OF STATE UPDATE METHODS                      !
  !--------------------------------------------------------------------------------!

  impure subroutine material_update(this, bootstrap)
    !! This subroutine updates the current estimate for the state of the material by numerically solving the diffusion equation.
    use bvp_m

    class(material),   intent(inout) :: this                       !! Material that will be updated
    logical, optional, intent(in)    :: bootstrap                  !! This flag prevents update posthooks
    type(bvp_sol)                    :: sol                        !  Workspace for the bvp_solver procedures
    type(propagator)                 :: a                          !  State at this energy at the left  interface
    type(propagator)                 :: b                          !  State at this energy at the right interface
    complex(wp)                      :: e                          !  Complex energy relative to the Thouless energy
    real(wp)                         :: u(32,size(this%location))  !  Representation of the retarded propagators
    real(wp)                         :: d(32,size(this%location))  !  Work array used to calculate the change in u(·,·)
    integer                          :: n                          !  Outer loop variable (current energy)
    integer                          :: m                          !  Inner loop variable (current location)

    ! Call the prehook method
    call this%update_prehook

    ! Only update the material state if it has a positive order
    if (this%order > 0) then
      ! Status information
      if (this%information >= 0) then
        block
          character(len=10) :: str
          write(str, '(2x,"[1m#",i0)') this%order
          write(stdout,'(a,a,20x)') str, trim(this%type_string)
        end block
      end if

      ! Reset the difference since last update to zero
      this%difference = 0.0_wp

      ! Loop over the discretized energy levels
      do n=size(this%energy),1,-1
        ! Status information
        if (this%information >= 0) then
          write(stdout,'(6x,a,1x,i4,1x,a,1x,i4,1x,a,f0.5,a1)',advance='no') &
            '[',n,'/',size(this%energy),']  E = ',this%energy(n), achar(13)
          flush(stdout)
        end if

        ! Convert all states at this energy level to real-valued state vectors
        do m=1,size(this%location)
          d(:,m) = this%propagator(n,m)
        end do

        ! If there is little difference between the previous state and this one, use the previous state as an initial guess
        if (this%difference < 0.05_wp) then
          u = d
        end if

        ! Calculate the complex energy (relative to the Thouless energy)
        e = cx(this%energy(n), this%scattering)/this%thouless

        ! Update the matrices used to evaluate boundary conditions
        if (associated(this%material_a)) then
          a = this%material_a%propagator(n,ubound(this%material_a%propagator,2))
        else
          a = propagator()
        end if

        if (associated(this%material_b)) then
          b = this%material_b%propagator(n,lbound(this%material_b%propagator,2))
        else
          b = propagator()
        end if

        ! Initialize bvp_solver
        sol = bvp_init(32, 16, this%location, u, max_num_subintervals=(size(this%location)*this%scaling))

        ! Solve the differential equation
        sol = bvp_solver(sol, ode, bc, method=this%method, error_control=this%control, &
                         tol=this%tolerance, trace=this%information, stop_on_fail=.false.)

        ! Check if the calculation succeeded
        if (sol%info /= 0) then
          call error('Failed to converge! This is usually because of an ill-posed problem, e.g. a system with no superconductors.')
        end if

        ! Use the results to update the state
        call bvp_eval(sol, this%location, u)
        do m=1,size(this%location)
          this%propagator(n,m) = u(:,m)
        end do

        ! Update the difference vector
        d = abs(u - d)

        ! Update the maximal difference since last update
        this%difference = max(this%difference,maxval(d))

        ! Clean up after bvp_solver
        call bvp_terminate(sol)
      end do

      ! Status information
      if (this%information >= 0) then
        write(stdout,'(6x,a,f10.8,a)') 'Max change: ',this%difference,'                                        '
        flush(stdout)
      end if
    end if

    ! Stop here if bootstrapping
    if (present(bootstrap)) then
      if (bootstrap) then
        return
      end if
    end if

    ! Call the posthook method
    call this%update_posthook
  contains
    pure subroutine ode(z, u, f)
      !! Definition of the differential equation u'=f(z,u).
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
      call this%interface_equation_a(a, g1, gt1, dg1, dgt1, r1, rt1)
      call this%interface_equation_b(b, g2, gt2, dg2, dgt2, r2, rt2)

      ! Pack the results into state vectors
      bca(1: 8) = r1
      bca(9:16) = rt1
      bcb(1: 8) = r2
      bcb(9:16) = rt2
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
      case("order")
        call evaluate(val, this%order)
        if (this%order > 16) then
          call error("The order of the material must be be maximum 16.")
        else if (this%order < 0) then
          call error("The order of the material must be be minimum 0.")
        end if
      case default
        call warning("Unknown option '" // key // "' ignored.")
    end select
  end subroutine

  impure subroutine material_save(this)
    !! Save a backup of the current material state.
    class(material), intent(inout) :: this

    if (.not. allocated(this%backup)) then
      allocate(this%backup(size(this%propagator,1),size(this%propagator,2)))
    end if
    this%backup = this%propagator
  end subroutine

  impure subroutine material_load(this)
    !! Load a backup of a previous material state.
    class(material), intent(inout) :: this
    integer                        :: info

    ! Load the saved propagators
    if (allocated(this%backup)) then
      this%propagator = this%backup
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
