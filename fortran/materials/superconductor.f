!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This module defines the data type 'superconductor', which models the physical state of a superconductor. The type is
!> a member of class(conductor), and thus inherits the internal structure and generic methods defined in conductor_m.

module superconductor_m
  use :: stdio_m
  use :: condmat_m
  use :: conductor_m
  use :: ferromagnet_m
  private

  ! Type declaration
  type, public, extends(ferromagnet) :: superconductor
    ! These parameters control the physical characteristics of the material 
    complex(wp), allocatable :: gap_history(:,:)     !! Superconducting order parameter as a function of location (backup of previously calculated gaps on the location mesh)
    complex(wp), allocatable :: gap_function(:)      !! Superconducting order parameter as a function of location (relative to the zero-temperature gap of a bulk superconductor)
    real(wp),    allocatable :: gap_location(:)      !! Location array for the gap function (required because we interpolate the gap to a higher resolution than the propagators)
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
    procedure                :: gap                 => superconductor_gap                 !! Return the superconducting order parameter at a given position

    ! These methods define miscellaneous utility functions
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
    call this % ferromagnet % construct

    ! Initialize superconductivity
    this % correlation = 1

    ! Enable a higher-order solver
    this % method = 6

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
    call this % ferromagnet % initialize

    ! Disable status messages
    info = this % information
    if (this % information >= 0) then
      this % information = -1
    end if

    ! Silently initialize the gap
    call this % update_gap

    ! Reenable status messages
    this % information = info

    ! Reset the iteration counter
    this % iteration = 0
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
      gap  = (this % gap(z))/(this % thouless)
      gapt = conjg(gap)

      ! Calculate the second derivatives of the Riccati parameters (superclass terms)
      call this % ferromagnet % diffusion_equation(p, e, z)

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
    call this % ferromagnet % kinetic_equation(Gp, R, z)

    ! Lookup the superconducting order parameter
    gap  = (this % gap(z))/(this % thouless)
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
    call this % ferromagnet % update_prehook

    ! Modify the type string
    this % type_string = color_green // 'SUPERCONDUCTOR' // color_none
  end subroutine

  impure subroutine superconductor_update_posthook(this)
    !! Updates the superconducting order parameter based on the propagators of the system.
    class(superconductor), intent(inout) :: this

    ! Call the superclass posthook
    call this % ferromagnet % update_posthook

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

    associate( raw_loc  => this % location,        &
               raw_val  => this % correlation,     &
               raw_bak  => this % gap_history,     &
               gap_loc  => this % gap_location,    &
               gap_val  => this % gap_function,    &
               r        => this % progressive,     &
               m  => lbound(this % gap_history,2), &
               n  => ubound(this % gap_history,2)  )

      ! Calculate difference from previous gap
      diff = mean(abs(raw_val(:) - raw_bak(:,n)))

      ! Linear mixing of the old and new solutions
      raw_val = r*raw_val + (1-r)*raw_bak(:,n)

      ! Save the calculated gap as backup
      raw_bak(:,m:n-1) = raw_bak(:,m+1:n)
      raw_bak(:, n   ) = raw_val

      ! Interpolate the gap as a function of position to a higher resolution
      gap_val = interpolate(raw_loc, raw_val, gap_loc)
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
    !! where g(Δ)=f(Δ)-Δ=0, and this can be done more efficiently using e.g. Newtons 
    !! method than a straight-forward fixpoint-iteration. Using Newtons method, we get:
    !!   Δ_{n+1} = Δ_{n} - g(Δ_n)/g'(Δ_n)
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

    ! Stop here if it is not yet time to boost
    if (this % iteration > 0) then
      return
    end if

    ! Boost the convergence using Steffensen's method
    d1 = f(2) - f(1)
    d2 = f(3) - 2*f(2) + f(1)
    g  = f(1) - d1**2/d2

    ! Avoid boosts near transparent interfaces
    if (this % transparent_a) then
      g(lbound(g,1)) = this % gap(0.0_wp)
    end if
    if (this % transparent_b) then
      g(ubound(g,1)) = this % gap(1.0_wp)
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

  pure function superconductor_gap(this, location) result(gap)
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
        call this % ferromagnet % conf(key, val)
    end select
  end subroutine
end module
