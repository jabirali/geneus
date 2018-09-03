!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This module defines the data type 'material',  which models the state of a physical material for a discretized range of
!> positions and energies. This is an abstract type, meaning that it is not intended to be instantiated on its own, but is
!> intended as a base type for physical materials like conductors, superconductors, and ferromagnets. In other words, this
!> type defines the essential data structures and program structure, while the derived subtypes will define actual physics.

module material_m
  use :: stdio_m
  use :: condmat_m
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
    logical                                   :: phaselock             =  .false.              !! Lock the center-of-mass phase?

    ! Material properties that affect the nonequilibrium state
    logical                                   :: nonequilibrium        =  .false.              !! Equilibrium?
    logical                                   :: transverse            =  .false.              !! Transverse potential gradients?
    real(wp)                                  :: voltage               =  0.00_wp              !! Voltage (eV/Δ₀)
    real(wp)                                  :: temperature           =  0.01_wp              !! Temperature (T/Tc)
    real(wp)                                  :: spinvoltage           =  0.00_wp              !! Spin-voltage (eVs/Δ₀)
    real(wp)                                  :: spintemperature       =  0.00_wp              !! Spin-temperature (Ts/Tc)
    real(wp),                  dimension(1:3) :: spinaxis              =  [0,0,1]              !! Spin quantization axis

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
    real(wp),                     allocatable :: magnetization(:,:)                            !! Magnetization due to exchange and Zeeman effects
    complex(wp),                  allocatable :: correlation(:)                                !! Superconducting pair-correlations

    ! Hybrid structures are modeled as double-linked lists
    integer                                   :: order           =  1                          !! Simulation priority of this material
    class(material),                  pointer :: material_a      => null()                     !! Material to the left  (default: vacuum)
    class(material),                  pointer :: material_b      => null()                     !! Material to the right (default: vacuum)

    ! Control parameters for the numerical solvers
    integer                                   :: iteration       =  0                          !! Used to keep track of selfconsistent iteration cycles
    integer                                   :: selfconsistency =  2                          !! Selfconsistency scheme (0 = none, 1 = fixpoint, 2 = boost)
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
    procedure(manipulate),           deferred :: construct                                     !! Construct  object
    procedure(initialize),           deferred :: initialize                                    !! Initialize object
    procedure(manipulate),           deferred :: update_prehook                                !! Executed before update
    procedure(manipulate),           deferred :: update_posthook                               !! Executed after  update
    procedure                                 :: update           => material_update           !! Calculate propagators
    procedure                                 :: update_diffusion => diffusion_update          !! Calculate propagators (in equilibrium)
    procedure                                 :: update_kinetic   => kinetic_update            !! Calculate propagators (nonequilibrium)

    ! These methods define the physical equations used by the update methods
    procedure(diffusion_equation),   deferred :: diffusion_equation                            !! Diffusion equation
    procedure(diffusion_equation_a), deferred :: diffusion_equation_a                          !! Boundary condition (left)
    procedure(diffusion_equation_b), deferred :: diffusion_equation_b                          !! Boundary condition (right)

    procedure(kinetic_equation),     deferred :: kinetic_equation                              !! Kinetic equation
    procedure(kinetic_equation_a),   deferred :: kinetic_equation_a                            !! Boundary condition (left)
    procedure(kinetic_equation_b),   deferred :: kinetic_equation_b                            !! Boundary condition (right)

    ! These methods define miscellaneous utility functions
    procedure                                 :: conf            => material_conf              !! Configures material parameters
    procedure                                 :: save            => material_save              !! Saves the state of the material
    procedure                                 :: load            => material_load              !! Loads the state of the material
  end type

  ! Declare submodule procedures
  interface
    module subroutine kinetic_update(this)
      class(material), target, intent(inout) :: this
    end subroutine

    module subroutine diffusion_update(this)
      class(material), target, intent(inout) :: this
    end subroutine
  end interface

  ! Declare subclass procedures
  abstract interface
    subroutine initialize(this)
      !! This interface is used for the deferred procedure initialize.
      import material, wp

      class(material), intent(inout) :: this
    end subroutine

    subroutine manipulate(this)
      !! This interface is used for the deferred procedures construct, update_prehook, and update_posthook.
      import material

      class(material), intent(inout) :: this
    end subroutine

    pure subroutine diffusion_equation(this, p, e, z)
      !! This interface is used for the deferred procedure diffusion_equation.
      import material, propagator, wp

      class(material),  intent(in)    :: this
      type(propagator), intent(inout) :: p
      complex(wp),      intent(in)    :: e
      real(wp),         intent(in)    :: z
    end subroutine

    pure subroutine diffusion_equation_a(this, p, a, r, rt)
      !! This interface is used for the deferred procedure diffusion_equation_a.
      import material, propagator, spin, wp

      class(material),          intent(in)    :: this
      type(propagator),         intent(in)    :: p, a
      type(spin),               intent(inout) :: r, rt
    end subroutine

    pure subroutine diffusion_equation_b(this, p, b, r, rt)
      !! This interface is used for the deferred procedure diffusion_equation_b.
      import material, propagator, spin, wp

      class(material),          intent(in)    :: this
      type(propagator),         intent(in)    :: p, b
      type(spin),               intent(inout) :: r, rt
    end subroutine

    pure subroutine kinetic_equation(this, Gp, R, z)
      !! This interface is used for the deferred procedure kinetic_equation.
      import material, propagator, wp

      class(material),                 intent(in)    :: this
      type(propagator),                intent(in)    :: Gp
      complex(wp), dimension(0:7,0:7), intent(inout) :: R
      real(wp),                        intent(in)    :: z
    end subroutine

    pure subroutine kinetic_equation_a(this, Gp, Ga, Cp, Ca)
      !! This interface is used for the deferred procedure kinetic_equation_a.
      import material, propagator, wp

      class(material),                 intent(in)  :: this
      type(propagator),                intent(in)  :: Gp, Ga
      complex(wp), dimension(0:7,0:7), intent(out) :: Cp, Ca
    end subroutine

    pure subroutine kinetic_equation_b(this, Gp, Gb, Cp, Cb)
      !! This interface is used for the deferred procedure kinetic_equation_b.
      import material, propagator, wp

      class(material),                 intent(in)  :: this
      type(propagator),                intent(in)  :: Gp, Gb
      complex(wp), dimension(0:7,0:7), intent(out) :: Cp, Cb
    end subroutine
  end interface
contains

  !--------------------------------------------------------------------------------!
  !                    IMPLEMENTATION OF STATE UPDATE METHODS                      !
  !--------------------------------------------------------------------------------!

  impure subroutine material_update(this, bootstrap)
    !! This subroutine updates the state of the material by solving the diffusion 
    !! equations for the equilibrium propagators, the kinetic equations for the 
    !! nonequilibrium propagators, and then calculating physical observables.
    class(material),   intent(inout) :: this      !! Material that will be updated
    logical, optional, intent(in)    :: bootstrap !! Disable calculation of observables

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
    if (this % nonequilibrium) then
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
      case("voltage")
        call evaluate(val, this%voltage)
      case("spinvoltage")
        call evaluate(val, this%spinvoltage)
      case("spintemperature")
        call evaluate(val, this%spintemperature)
      case("spinaxis")
        call evaluate(val, this%spinaxis)
        this%spinaxis = unitvector(this%spinaxis)
      case("transverse")
        call evaluate(val, this%transverse)
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
      case("iteration")
        call evaluate(val, this%iteration)
      case("selfconsistency")
        call evaluate(val, this%selfconsistency)
        if (this % selfconsistency < 0 .or. this % selfconsistency > 2) then
          call error("The selfconsistency scheme should be in the range [0,2].")
        end if
      case("phaselock")
        call evaluate(val, this%phaselock)
      case("equilibrium")
        call evaluate(val, this%nonequilibrium)
        this % nonequilibrium = .not. this % nonequilibrium
      case("nonequilibrium")
        call evaluate(val, this%nonequilibrium)
      case default
        call error("Unknown material option '" // key // "'.")
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

    ! Load saved propagators
    call this % propagator % load(this % backup)

    ! Disable status messages
    info = this%information
    if (this%information >= 0) then
      this%information = -1
    end if

    ! Silently recalculate derived properties
    call this % update_posthook

    ! Reenable status messages
    this % information = info
  end subroutine
end module
