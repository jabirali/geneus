!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This module defines the data type 'superconductor', which models the physical state of a superconductor. The type is
!> a member of class(conductor), and thus inherits the internal structure and generic methods defined in conductor_m.

module superconductor_m
  use :: stdio_m
  use :: condmat_m
  use :: conductor_m
  private

  ! Type declaration
  type, public, extends(conductor) :: superconductor
    ! These parameters control the physical characteristics of the material 
    complex(wp), allocatable :: gap_history(:,:)     !! Superconducting order parameter as a function of location (backup of previously calculated gaps on the location mesh)
    complex(wp), allocatable :: gap_function(:)      !! Superconducting order parameter as a function of location (relative to the zero-temperature gap of a bulk superconductor)
    real(wp),    allocatable :: gap_location(:)      !! Location array for the gap function (required because we interpolate the gap to a higher resolution than the propagators)
    integer                  :: iteration            !! Used to count where in the selfconsistent iteration cycle we are
  contains
    ! These methods define the class(material) interface
    procedure                :: construct          => superconductor_construct            !! Construct  propagators
    procedure                :: initialize         => superconductor_initialize           !! Initialize propagators

    ! These methods contain the equations that describe superconductors
    procedure                :: diffusion_equation  => superconductor_diffusion_equation  !! Diffusion equation
    procedure                :: kinetic_equation    => superconductor_kinetic_equation    !! Kinetic equation
    procedure                :: update_gap          => superconductor_update_gap          !! Calculate the superconducting order parameter
    procedure                :: update_boost        => superconductor_update_boost        !! Boost the convergence of the order parameter (Steffensen's method)
    procedure                :: update_prehook      => superconductor_update_prehook      !! Update the internal variables before calculating the propagators
    procedure                :: update_posthook     => superconductor_update_posthook     !! Update the superconducting order parameter from  the propagators

    ! These methods are used to access and mutate the parameters
    procedure                :: set_gap             => superconductor_set_gap             !! Update the superconducting order parameter from a given scalar
    procedure                :: get_gap             => superconductor_get_gap             !! Return the superconducting order parameter at a given position

    ! These methods define miscellaneous utility functions
    procedure                :: load                => superconductor_load                !! Load the state of a superconductor
    procedure                :: conf                => superconductor_conf                !! Configure material parameters
  end type
contains

  !--------------------------------------------------------------------------------!
  !                        IMPLEMENTATION OF CONSTRUCTORS                          !
  !--------------------------------------------------------------------------------!

  impure subroutine superconductor_construct(this)
    !! Constructs a superconducting material that is initialized to a superconducting state.
    class(superconductor), intent(inout) :: this

    ! Call the superclass constructor
    call this % conductor % construct

    ! Initialize superconductivity
    this % correlation = 1

    ! Allocate interpolation functions
    allocate(this % gap_history(size(this%location),1:3))
    allocate(this % gap_location(4096 * size(this%location)))
    allocate(this % gap_function(size(this%gap_location)))
    call linspace(this % gap_location, this % location(1), this % location(size(this % location)))
  end subroutine

  impure subroutine superconductor_initialize(this)
    !! Redefine the default initializer.
    class(superconductor), intent(inout) :: this
    integer                              :: info

    ! Call the superclass initializer
    call this % conductor % initialize()

    ! Disable status messages
    info = this%information
    if (this%information >= 0) then
      this%information = -1
    end if

    ! Silently initialize the gap
    call this % update_gap

    ! Reenable status messages
    this%information = info
  end subroutine

  !--------------------------------------------------------------------------------!
  !                   IMPLEMENTATION OF SUPERCONDUCTOR METHODS                     !
  !--------------------------------------------------------------------------------!

  pure subroutine superconductor_diffusion_equation(this, p, e, z)
    !! Use the diffusion equation to calculate the second derivatives of the Riccati parameters at point z.
    class(superconductor), intent(in)    :: this
    complex(wp),           intent(in)    :: e
    real(wp),              intent(in)    :: z
    type(propagator),      intent(inout) :: p
    complex(wp)                          :: gap, gapt

    associate(  N => p % N,     Nt => p % Nt,  &
                g => p % g,     gt => p % gt,  &
               dg => p % dg,   dgt => p % dgt, &
              d2g => p % d2g, d2gt => p % d2gt )

      ! Lookup the superconducting order parameter
      gap  = this%get_gap(z)/this%thouless
      gapt = conjg(gap)

      ! Calculate the second derivatives of the Riccati parameters (conductor terms)
      call this % conductor % diffusion_equation(p, e, z)

      ! Calculate the second derivatives of the Riccati parameters (superconductor terms)
      d2g  = d2g  - gap  * pauli2 + gapt * g  * pauli2 * g
      d2gt = d2gt + gapt * pauli2 - gap  * gt * pauli2 * gt
    end associate
  end subroutine

  pure subroutine superconductor_kinetic_equation(this, Gp, R, z)
    !! Calculate the self-energies in the kinetic equation.
    class(superconductor),           intent(in)    :: this
    type(propagator),                intent(in)    :: Gp
    complex(wp), dimension(0:7,0:7), intent(inout) :: R
    real(wp),                        intent(in)    :: z
    complex(wp)                                    :: gap, gapt
    type(nambu)                                    :: S

    ! Call the superclass kinetic equation
    call this % conductor % kinetic_equation(Gp, R, z)

    ! Lookup the superconducting order parameter
    gap  = (this % get_gap(z))/(this % thouless)
    gapt = conjg(gap)

    ! Construct the self-energy matrix
    S % matrix(1,4) = +gap
    S % matrix(2,3) = -gap
    S % matrix(3,2) = +gapt
    S % matrix(4,1) = -gapt

    ! Calculate the self-energy contribution
    R = R + Gp % selfenergy1(S)
  end subroutine

  impure subroutine superconductor_update_prehook(this)
    !! Code to execute before running the update method of a class(superconductor) object.
    class(superconductor), intent(inout) :: this

    ! Call the superclass prehook
    call this % conductor % update_prehook

    ! Modify the type string
    this % type_string = color_green // 'SUPERCONDUCTOR' // color_none
  end subroutine

  impure subroutine superconductor_update_posthook(this)
    !! Updates the superconducting order parameter based on the propagators of the system.
    class(superconductor), intent(inout) :: this

    ! Call the superclass posthook
    call this % conductor % update_posthook

    ! Update the superconducting gap using fixpoint-iteration
    if (this % selfconsistency >= 1) then
      call this % update_gap
    end if

    ! Boost the superconducting gap using Steffensen's method
    if (this % selfconsistency >= 2) then
      call this % update_boost
    end if
  end subroutine

  impure subroutine superconductor_update_gap(this)
    !! Interpolate the superconducting correlations Δ(z) to a higher resolution,
    !! to make the calculations more stable near strong ferromagnetic materials.
    class(superconductor), intent(inout)        :: this       !! Superconductor object
    real(wp)                                    :: diff       !! Change between iterations

    ! Interpolate the gap as a function of position to a higher resolution
    this % gap_function = interpolate(this % location, this % correlation, this % gap_location)

    ! Save the calculated gap as backup
    associate( b => this % gap_history, m => lbound(this % gap_history,2), n => ubound(this % gap_history,2) )
      b(:,m:n-1) = b(:,m+1:n)
      b(:,  n  ) = this % correlation
      diff       = mean(abs(b(:,n) - b(:,n-1)))
    end associate

    ! Status information
    if (this%information >= 0 .and. this%order > 0) then
      write(stdout,'(6x,a,f10.8,a,10x)') 'Gap change: ', diff
      flush(stdout)
    end if
  end subroutine

  impure subroutine superconductor_update_boost(this)
    !! Boost the convergence of the order parameter using Steffensen's method.
    !!
    !! The basic idea is that a selfconsistent solution of the Usadel equations can be
    !! regarded as a fixpoint iteration problem Δ=f(Δ), where the function f consists
    !! of solving the diffusion and gap equations one time. We're seeking the point
    !! where f'(Δ)=0, and this can be done more efficiently using e.g. Newtons method
    !! than a straight-forward fixpoint-iteration. Using Newtons method, we get:
    !!   Δ_{n+1} = Δ_{n} - f'(Δ_n)/f''(Δ_n)
    !! Using a finite-difference approximation for the derivatives, we arrive at the
    !! Steffensen iteration scheme, which yields an improved 2nd-order convergence:
    !!   Δ_{n+3} = Δ_{n} - (Δ_{n+1} - Δ_{n})²/(Δ_{n+2} - 2Δ_{n+1} + Δ_{n})
    !! In most of my tests, the simulation time is then reduced by a factor 2x to 5x.
    !!
    !! @NOTE:
    !!   I have also experimented with several higher-order methods, including methods
    !!   that utilize the 3rd derivative or perform multiple successive boosts. However, 
    !!   my experience so far is that these are less stable than Steffensen's method,
    !!   and often converged more slowly. These methods have therefore been discarded.

    class(superconductor), intent(inout)        :: this
    complex(wp), dimension(size(this%location)) :: g, d1, d2
    logical                                     :: u

    ! Update the iterator
    this % iteration = modulo(this % iteration + 1, 8)

    ! Right after a convergence boost, a 6th-order Runge-Kutta algorithm is the 
    ! most efficient alternative for solving the diffusion equations. Otherwise,
    ! a 4th-order algorithm is sufficient, and spends less time per iteration.
    if (this % iteration == 0) then
      this % method = 6
    else
      this % method = 4
    end if

    ! Stop here if it is not yet time to boost
    if (this % iteration > 0) then
      return
    end if

    ! Boost the convergence using Steffensen's method
    d1 = f(2) - f(1)
    d2 = f(3) - 2*f(2) + f(1)
    g  = f(1) - d1**2/d2

    ! Abort now if the boost is numerically unstable
    if (any(abs(d2) < 1e-10)) then
      return
    end if

    ! Interpolate the gap as a function of position to a higher resolution
    this % gap_function = interpolate(this % location, g, this % gap_location)

    ! Perform one extra update if necessary
    u = .false.
    if (associated(this % material_a)) then
      u = u .or. (this % material_a % order > 0)
    end if
    if (associated(this % material_b)) then
      u = u .or. (this % material_b % order > 0)
    end if
    if (u) then
      call this % update(bootstrap = .true.)
    end if

    ! Status information
    if (this%information >= 0 .and. this%order > 0) then
      write(stdout,'(6x,a,f10.8,a,10x)') 'Gap boost:  ', mean(abs(g-f(3)))
      flush(stdout)
    end if
  contains
    pure function f(n) result(r)
      ! Define an accessor for the gap after n iterations.
      integer, intent(in)                              :: n
      complex(wp), dimension(size(this%gap_history,1)) :: r

      r = this % gap_history(:,n)
    end function
  end subroutine

  !--------------------------------------------------------------------------------!
  !                    IMPLEMENTATION OF GETTERS AND SETTERS                       !
  !--------------------------------------------------------------------------------!

  pure subroutine superconductor_set_gap(this, gap)
    !! Updates the superconducting order parameter from a scalar.
    class(superconductor), intent(inout) :: this
    complex(wp),           intent(in)    :: gap

    this % correlation  = gap
    this % gap_function = gap
    this % gap_history  = gap
    this % iteration    = 0
  end subroutine

  pure function superconductor_get_gap(this, location) result(gap)
    !! Returns the superconducting order parameter at the given location.
    class(superconductor), intent(in) :: this
    real(wp),              intent(in) :: location
    complex(wp)                       :: gap
    integer                           :: n, m

    associate(f => this%gap_function, fp => gap, p => location)
      ! Calculate the index corresponding to the given location
      m = size(f) - 1
      n = floor(p*m + 1)

      ! Interpolate the superconducting order parameter at that point
      if (n <= 1) then
        fp = f(1)
      else
        fp = f(n-1) + (f(n)-f(n-1))*(p*m-(n-2))
      end if
    end associate
  end function

  !--------------------------------------------------------------------------------!
  !                      IMPLEMENTATION OF UTILITY METHODS                         !
  !--------------------------------------------------------------------------------!

  impure subroutine superconductor_conf(this, key, val)
    !! Configure a material property based on a key-value pair.
    use :: evaluate_m

    class(superconductor), intent(inout) :: this
    character(*),          intent(in)    :: key
    character(*),          intent(in)    :: val

    select case(key)
      case default
        call this % conductor % conf(key, val)
    end select
  end subroutine

  impure subroutine superconductor_load(this)
    !! Load a backup of a previous material state.
    use :: material_m
    class(superconductor), intent(inout) :: this

    ! Call the superclass prehook
    call material_load(this)

    ! Reset the iteration counter
    this % iteration = 0
  end subroutine
end module
