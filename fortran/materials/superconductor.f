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
    complex(wp), allocatable :: gap_backup(:,:)    ! Superconducting order parameter as a function of location (backup of previously calculated gaps on the location mesh)
    complex(wp), allocatable :: gap_function(:)    ! Superconducting order parameter as a function of location (relative to the zero-temperature gap of a bulk superconductor)
    real(wp),    allocatable :: gap_location(:)    ! Location array for the gap function (required because we interpolate the gap to a higher resolution than the propagators)
    real(wp)                 :: coupling  =  0     ! BCS coupling constant that defines the strength of the superconductor (dimensionless)
    integer                  :: boost     =  3     ! What kind of selfconsistent iteration scheme to use (0 = fixpoint, 1 = Steffensen's method, 2 = Ostrowski's method, 3 = Cordero's method)
    integer                  :: iteration = -1     ! Used to count where in the selfconsistent iteration cycle we are
  contains
    ! These methods contain the equations that describe superconductors
    procedure                :: init                => superconductor_init                ! Initializes the propagators
    procedure                :: diffusion_equation  => superconductor_diffusion_equation  ! Defines the Usadel diffusion equation
    procedure                :: update_gap          => superconductor_update_gap          ! Calculates the superconducting order parameter
    procedure                :: update_boost        => superconductor_update_boost        ! Boosts the convergence of the order parameter
    procedure                :: update_prehook      => superconductor_update_prehook      ! Updates the internal variables before calculating the propagators
    procedure                :: update_posthook     => superconductor_update_posthook     ! Updates the superconducting order parameter from  the propagators

    ! These methods are used to access and mutate the parameters
    procedure                :: set_gap             => superconductor_set_gap             ! Updates the superconducting order parameter from a given scalar
    procedure                :: get_gap             => superconductor_get_gap             ! Returns the superconducting order parameter at a given position
    procedure                :: get_gap_mean        => superconductor_get_gap_mean        ! Returns the superconducting order parameter averaged over the material

    ! These methods define miscellaneous utility functions
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
    allocate(this%gap_backup(size(this%location),0:6))
    allocate(this%gap_location(4096 * size(this%location)))
    allocate(this%gap_function(size(this%gap_location)))
    call linspace(this%gap_location, this%location(1), this%location(size(this%location)))
    call this%set_gap( (1.0_wp,0.0_wp) )

    ! Initialize the coupling constant
    if (this%energy(size(this%energy)) > 0) then
      this%coupling = 1/acosh(this%energy(size(this%energy)))
    end if
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

    ! Update the superconducting gap
    call this%update_gap

    ! Boost the convergence if possible
    call this%update_boost
  end subroutine

  impure subroutine superconductor_update_gap(this)
    ! Update the superconducting gap if the BCS coupling constant is nonzero.
    use :: calculus_m

    class(superconductor), intent(inout) :: this         ! Superconductor object that will be updated
    complex(wp),           allocatable   :: gap_e(:)     ! Used to calculate the superconducting order parameter (as a function of energy)
    complex(wp),           allocatable   :: gap_z(:)     ! Used to calculate the superconducting order parameter (as a function of position)
    complex(wp)                          :: diff         ! Change in the superconducting order parameter mean
    complex(wp)                          :: singlet      ! Singlet component of the anomalous propagators
    integer                              :: n, m         ! Loop variables

    if (abs(this%coupling) > eps) then
      ! Allocate workspace memory
      allocate(gap_e(size(this%energy)))
      allocate(gap_z(size(this%location)))

      ! Calculate the mean superconducting order parameter
      diff = this%get_gap_mean()

      ! Iterate over the stored propagators
      do n = 1,size(this%location)
        do m = 1,size(this%energy)
          ! Calculate the singlet components of the anomalous propagators
          singlet = ( this%propagator(m,n)%singlet() - conjg(this%propagator(m,n)%singlett()) )/2.0_wp

          ! Calculate the gap equation integrand and store it in an array
          gap_e(m)  = singlet * this%coupling * tanh(0.8819384944310228_wp * this%energy(m)/this%temperature)
        end do

        ! Interpolate and integrate the results, and update the superconducting order parameter
        gap_z(n) = integrate(this%energy, gap_e, 1e-6_wp, this%energy(ubound(this%energy,1)))
      end do

      ! Save the calculated gap as backup
      associate( b => this % gap_backup, m => lbound(b,2), n => ubound(b,2) )
        b(:,m:n-1) = b(:,m+1:n)
        b(:,  n  ) = gap_z
      end associate

      ! Interpolate the gap as a function of position to a higher resolution
      this % gap_function = interpolate(this % location, gap_z, this % gap_location)

      ! Calculate the difference in mean superconducting order parameter
      diff = this%get_gap_mean() - diff

      ! Status information
      if (this%information >= 0 .and. this%order > 0) then
        write(stdout,'(6x,a,f10.8,a,10x)') 'Gap change: ',abs(diff)
        flush(stdout)
      end if

      ! Deallocate workspace memory
      deallocate(gap_e)
      deallocate(gap_z)
    end if
  end subroutine

  impure subroutine superconductor_update_boost(this)
    !! Boost the convergence of the order parameter using either Steffensen's method
    !! (2nd-order), Ostrowski's method (4th-order), or Cordero's method (8th-order).
    !!
    !! The basic idea is that a selfconsistent solution of the Usadel equations can be
    !! regarded as a fixpoint iteration problem Δ=f(Δ), where the function f consists
    !! of solving the diffusion and gap equations one time. We're seeking the point
    !! where f'(Δ)=0, and this can be done more efficiently using e.g. Newtons method
    !! than a straight-forward fixpoint-iteration. Using Newtons method, we get:
    !!   Δ_{n+1} = Δ_{n} - f'(Δ_n)/f''(Δ_n)
    !! Using a finite-difference approximation for the derivatives, we arrive at the
    !! Steffensen iteration scheme, which yields an improved 2nd-order convergence:
    !!   Δ_{n+1} = Δ_{n-2} - (Δ_{n-1} - Δ_{n-2})/(Δ_{n} - 2Δ_{n-1} + Δ_{n-2})
    !! Ostrowski's and Cordero's methods further improve the convergence by adding 
    !! more steps to the iteration scheme. See e.g. "A stable family with high order
    !! of convergence for solving nonlinear equations" by Cordero et al. (2015).

    use :: calculus_m

    class(superconductor), intent(inout)        :: this
    complex(wp), dimension(size(this%location)) :: g

    associate(g0 => this % gap_backup(:,0), &
              g1 => this % gap_backup(:,1), &
              g2 => this % gap_backup(:,2), &
              g3 => this % gap_backup(:,3), &
              g4 => this % gap_backup(:,4), &
              g5 => this % gap_backup(:,5), &
              g6 => this % gap_backup(:,6)  )

      ! Update the iterator
      this % iteration = this % iteration + 1

      ! Control the boost pattern
      select case (this % iteration)
        case (:5)
          ! Switch to a 4th-order Runge-Kutta method
          this % method = 4

          ! Reset the iterator
          if (this % boost <= 0) then
            this % iteration = 0
          end if

          ! Perform a regular fixpoint iteration
          return
        case (6)
          ! Switch to a 6th-order Runge-Kutta method
          this % method = 6

          ! Reset the iterator
          if (this % boost == 1) then
            this % iteration = 0
          end if

          ! Perform the 1st boost with Steffensen's method (2nd order)
          g = x(4) - d1(4)/d2(4,5)
        case (7)
          ! Perform the 2nd boost with Ostrowski's method (4th order)
          g = x(2) - (1 + d1(5)/(d1(2)-2*d1(5))) * d1(2)/d2(2,3)

          ! Reset the iterator
          if (this % boost == 2) then
            this % iteration = 0
          end if
        case (8)
          ! Perform the 3rd boost with Cordero's method (8th order)
          g = x(5) - d1(5)/(d2(5,3) + 2*(x(5)-x(3))*d3(5,3,0) - (x(5)-x(3))*d3(3,0,0))

          ! Reset the iterator
          this % iteration = 0
      end select

      ! Interpolate the gap as a function of position to a higher resolution
      this % gap_function = interpolate(this % location, g, this % gap_location)

      ! Perform one extra update
      call this % update(bootstrap = .true.)

      ! Status information
      if (this%information >= 0 .and. this%order > 0) then
        write(stdout,'(6x,a,f10.8,a,10x)') 'Gap boost:  ',sum(abs(g - g6))/size(g)
        flush(stdout)
      end if

      ! Make backups of the gaps
      g0 = g1
      g1 = g2
      g2 = g3
      g3 = g4
      g4 = g5
      g5 = g6
      g6 = g
    end associate
  contains
    pure function x(n) result(r)
      ! Define an accessor method for the gap value after n iterations.
      integer, intent(in)                             :: n
      complex(wp), dimension(size(this%gap_backup,1)) :: r

      r = this % gap_backup(:,n)
    end function

    pure function d1(n) result(r)
      ! Estimate the finite-difference 1st-derivative between iterations n and n+1.
      integer, intent(in)                             :: n
      complex(wp), dimension(size(this%gap_backup,1)) :: r

      r = x(n+1) - x(n)
    end function

    pure function d2(n,m) result(r)
      ! Estimate the finite-difference 2nd-derivative between iterations n and m.
      integer, intent(in)                             :: n, m
      complex(wp), dimension(size(this%gap_backup,1)) :: r

      if (n /= m) then
        r = (d1(m)   - d1(n))/(x(m)   - x(n))
      else
        r = (d1(n+1) - d1(n))/(x(n+1) - x(n))
      end if
    end function

    pure function d3(n,m,k) result(r)
      ! Estimate the finite-difference 3rd-derivative between iterations n, m, k.
      integer, intent(in)                             :: n, m, k
      complex(wp), dimension(size(this%gap_backup,1)) :: r

      r = (d2(m,k) - d2(n,m))/(x(k) - x(n))
    end function
  end subroutine

  !--------------------------------------------------------------------------------!
  !                    IMPLEMENTATION OF GETTERS AND SETTERS                       !
  !--------------------------------------------------------------------------------!

  pure subroutine superconductor_set_gap(this, gap)
    ! Updates the superconducting order parameter from a scalar.
    class(superconductor), intent(inout) :: this
    complex(wp),           intent(in)    :: gap

    this%gap_function = gap
    this%gap_backup   = gap
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

  pure function superconductor_get_gap_mean(this) result(gap)
    ! Returns the superconducting order parameter average in the material.
    class(superconductor), intent(in)  :: this
    complex(wp)                        :: gap

    gap = sum(this%gap_function)/max(1,size(this%gap_function)) 
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
    real(wp)                             :: tmp

    select case(key)
      case("boost")
        call evaluate(val, this%boost)
        if (this%boost < 0 .or. this%boost > 3) then
          call error("The 'boost' parameter must be in the range [0,3].")
        end if
      case("coupling")
        call evaluate(val, this%coupling)
      case ('gap')
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
end module
