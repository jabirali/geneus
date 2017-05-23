!> Author:   Jabir Ali Ouassou
!> Date:     2015-07-17
!> Category: Materials
!>
!> This module defines the data type 'superconductor', which models the physical state of a superconductor. The type is
!> a member of class(conductor), and thus inherits the internal structure and generic methods defined in conductor_m.

module superconductor_m
  use :: stdio_m
  use :: math_m
  use :: spin_m
  use :: conductor_m
  private

  ! Type declaration
  type, public, extends(conductor) :: superconductor
    ! These parameters control the physical characteristics of the material 
    complex(wp), allocatable :: gap_history(:,:)     ! Superconducting order parameter as a function of location (backup of previously calculated gaps on the location mesh)
    complex(wp), allocatable :: gap_function(:)      ! Superconducting order parameter as a function of location (relative to the zero-temperature gap of a bulk superconductor)
    real(wp),    allocatable :: gap_location(:)      ! Location array for the gap function (required because we interpolate the gap to a higher resolution than the propagators)
    integer                  :: selfconsistency = 2  ! What kind of selfconsistency scheme to use (0 = none, 1 = fixpoint-iteration, 2 = Steffensen's method)
    integer                  :: iteration            ! Used to count where in the selfconsistent iteration cycle we are
  contains
    ! These methods contain the equations that describe superconductors
    procedure                :: init                => superconductor_init                ! Initializes the propagators
    procedure                :: diffusion_equation  => superconductor_diffusion_equation  ! Defines the Usadel diffusion equation
    procedure                :: update_gap          => superconductor_update_gap          ! Calculates the superconducting order parameter
    procedure                :: update_boost        => superconductor_update_boost        ! Boosts the convergence of the order parameter (Steffensen's method)
    procedure                :: update_prehook      => superconductor_update_prehook      ! Updates the internal variables before calculating the propagators
    procedure                :: update_posthook     => superconductor_update_posthook     ! Updates the superconducting order parameter from  the propagators

    ! These methods are used to access and mutate the parameters
    procedure                :: set_gap             => superconductor_set_gap             ! Updates the superconducting order parameter from a given scalar
    procedure                :: get_gap             => superconductor_get_gap             ! Returns the superconducting order parameter at a given position

    ! These methods define miscellaneous utility functions
    procedure                :: load                => superconductor_load                ! Loads the state of a superconductor
    procedure                :: conf                => superconductor_conf                ! Configures material parameters
  end type

  ! Type constructor
  interface superconductor
    module procedure superconductor_construct
  end interface
contains

  !--------------------------------------------------------------------------------!
  !                        IMPLEMENTATION OF CONSTRUCTORS                          !
  !--------------------------------------------------------------------------------!

  function superconductor_construct() result(this)
    ! Constructs a superconducting material that is initialized to a superconducting state.
    use :: calculus_m
    type(superconductor) :: this

    ! Call the superclass constructor
    this%conductor = conductor()

    ! Initialize the order parameter
    allocate(this%gap_history(size(this%location),1:3))
    allocate(this%gap_location(4096 * size(this%location)))
    allocate(this%gap_function(size(this%gap_location)))
    call linspace(this%gap_location, this%location(1), this%location(size(this%location)))
    call this%init( (1.0_wp,0.0_wp) )
  end function

  pure subroutine superconductor_init(this, gap)
    ! Redefine the default initializer.
    class(superconductor), intent(inout) :: this
    complex(wp),           intent(in)    :: gap

    ! Call the superclass initializer
    call this%conductor%init(gap)

    ! Update the superconducting gap
    call this%set_gap(gap)
  end subroutine

  !--------------------------------------------------------------------------------!
  !                   IMPLEMENTATION OF SUPERCONDUCTOR METHODS                     !
  !--------------------------------------------------------------------------------!

  pure subroutine superconductor_diffusion_equation(this, e, z, g, gt, dg, dgt, d2g, d2gt)
    ! Use the diffusion equation to calculate the second derivatives of the Riccati parameters at point z.
    class(superconductor), intent(in)    :: this
    complex(wp),           intent(in)    :: e
    real(wp),              intent(in)    :: z
    type(spin),            intent(in)    :: g, gt, dg, dgt
    type(spin),            intent(inout) :: d2g, d2gt
    complex(wp)                          :: gap, gapt

    ! Lookup the superconducting order parameter
    gap  = this%get_gap(z)/this%thouless
    gapt = conjg(gap)

    ! Calculate the second derivatives of the Riccati parameters (conductor terms)
    call this%conductor%diffusion_equation(e, z, g, gt, dg, dgt, d2g, d2gt)

    ! Calculate the second derivatives of the Riccati parameters (superconductor terms)
    d2g  = d2g  - gap  * pauli2 + gapt * g  * pauli2 * g
    d2gt = d2gt + gapt * pauli2 - gap  * gt * pauli2 * gt
  end subroutine

  impure subroutine superconductor_update_prehook(this)
    ! Code to execute before running the update method of a class(superconductor) object.
    class(superconductor), intent(inout) :: this

    ! Call the superclass prehook
    call this%conductor%update_prehook

    ! Modify the type string
    this%type_string = color_green // 'SUPERCONDUCTOR' // color_none
    if (allocated(this%spinorbit))       this%type_string = trim(this%type_string) // color_cyan   // ' [SOC]' // color_none
    if (allocated(this%magnetization_a)) this%type_string = trim(this%type_string) // color_purple // ' [SAL]' // color_none
    if (allocated(this%magnetization_b)) this%type_string = trim(this%type_string) // color_purple // ' [SAR]' // color_none
  end subroutine

  impure subroutine superconductor_update_posthook(this)
    ! Updates the superconducting order parameter based on the propagators of the system.
    class(superconductor), intent(inout) :: this

    ! Call the superclass posthook
    call this%conductor%update_posthook

    ! Update the superconducting gap using fixpoint-iteration
    if (this % selfconsistency >= 1) then
      call this%update_gap
    end if

    ! Boost the superconducting gap using Steffensen's method
    if (this % selfconsistency >= 2) then
      call this%update_boost
    end if
  end subroutine

  impure subroutine superconductor_update_gap(this)
    !! Calculate the superconducting gap Δ(z) from the propagators using a selfconsistency equation.
    !! @TODO: The tanh(...) has to be generalized for future nonequilibrium calculations.
    use :: calculus_m

    class(superconductor), intent(inout)        :: this       ! Superconductor object that will be updated
    complex(wp), dimension(size(this%energy))   :: gap_e      ! Used to calculate the order parameter (as a function of energy)
    complex(wp), dimension(size(this%location)) :: gap_z      ! Used to calculate the order parameter (as a function of position)
    complex(wp), dimension(0:3)                 :: f, ft      ! Singlet/triplet decomposition of the anomalous propagators
    real(wp)                                    :: coupling   ! The BCS coupling that gives rise to superconductivity
    real(wp)                                    :: diff       ! Change of the gap between two iterations
    integer                                     :: n, m       ! Loop variables

    ! Calculate the appropriate coupling from the cutoff energy
    coupling = 1/acosh(this%energy(size(this%energy)))

    ! Iterate over the stored propagators
    do n = 1,size(this%location)
      do m = 1,size(this%energy)
        ! Perform a singlet/triplet decomposition of the anomalous propagators
        call this % propagator(m,n) % decompose(f = f, ft = ft)

        ! Calculate the gap equation integrand and store it in an array
        gap_e(m)  = ((f(0) - conjg(ft(0)))/2) * coupling * tanh(0.8819384944310228_wp * this%energy(m)/this%temperature)
      end do

      ! Interpolate and integrate the results, and update the superconducting order parameter
      gap_z(n) = integrate(this%energy, gap_e, 1e-6_wp, this%energy(ubound(this%energy,1)))
    end do

    ! Interpolate the gap as a function of position to a higher resolution
    this % gap_function = interpolate(this % location, gap_z, this % gap_location)

    ! Save the calculated gap as backup
    associate( b => this % gap_history, m => lbound(this % gap_history,2), n => ubound(this % gap_history,2) )
      b(:,m:n-1) = b(:,m+1:n)
      b(:,  n  ) = gap_z
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
    !!   Δ_{n+3} = Δ_{n} - (Δ_{n+1} - Δ_{n})/(Δ_{n+2} - 2Δ_{n+1} + Δ_{n})
    !! In most of my tests, the simulation time is then reduced by a factor 2x to 5x.
    !!
    !! @NOTE:
    !!   I have also experimented with several higher-order methods, including methods
    !!   that utilize the 3rd derivative or perform multiple successive boosts. However, 
    !!   my experience so far is that these are less stable than Steffensen's method,
    !!   and often converged more slowly. These methods have therefore been discarded.

    use :: calculus_m

    class(superconductor), intent(inout)        :: this
    complex(wp), dimension(size(this%location)) :: g
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

    ! Boost the convergence using Steffensen's method when the iterator resets
    if (this % iteration == 0) then
      g = f(1) - (f(2)-f(1))**2/(f(3)-2*f(2)+f(1))
    else
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
    ! Updates the superconducting order parameter from a scalar.
    class(superconductor), intent(inout) :: this
    complex(wp),           intent(in)    :: gap

    this % gap_function = gap
    this % gap_history  = gap
    this % iteration    = 0
  end subroutine

  pure function superconductor_get_gap(this, location) result(gap)
    ! Returns the superconducting order parameter at the given location.
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
      case("selfconsistency")
        call evaluate(val, this%selfconsistency)
        if (this % selfconsistency < 0 .or. this % selfconsistency > 2) then
          call error("The selfconsistency scheme should be in the range [0,2].")
        end if
      case ("gap")
        block
          integer  :: index
          real(wp) :: gap, phase
          index = scan(val,',')
          if (index > 0) then
            call evaluate(val(1:index-1), gap)
            call evaluate(val(index+1: ), phase)
          else
            call evaluate(val, gap)
            phase = 0
          end if
          call this % init( gap = gap*exp((0.0,1.0)*pi*phase) )
        end block
      case default
        call this%conductor%conf(key, val)
    end select
  end subroutine

  impure subroutine superconductor_load(this)
    ! Load a backup of a previous material state.
    use :: material_m
    class(superconductor), intent(inout) :: this

    ! Call the superclass prehook
    call material_load(this)

    ! Reset the iteration counter
    this % iteration = 0
  end subroutine
end module
