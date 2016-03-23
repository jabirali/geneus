! This module defines the data type 'superconductor', which models the physical state of a superconductor. The type is
! a member of class(conductor), and thus inherits the internal structure and generic methods defined in mod_conductor.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-17
! Updated: 2016-03-23

module mod_superconductor
  use mod_stdio
  use mod_math
  use mod_spin
  use mod_conductor
  implicit none
  private

  ! Type declaration
  type, public, extends(conductor) :: superconductor
    ! These parameters control the physical characteristics of the material 
    complex(wp), allocatable :: gap(:)             ! Superconducting order parameter as a function of position (relative to the zero-temperature gap of a bulk superconductor)
    real(wp)                 :: coupling = 0.00_wp ! BCS coupling constant that defines the strength of the superconductor (dimensionless)
  contains
    ! These methods contain the equations that describe superconductors
    procedure                :: init                => superconductor_init                ! Initializes the propagators
    procedure                :: diffusion_equation  => superconductor_diffusion_equation  ! Defines the Usadel diffusion equation
    procedure                :: update_gap          => superconductor_update_gap          ! Calculates the superconducting order parameter
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

  pure function superconductor_construct() result(this)
    ! Constructs a superconducting material that is initialized to a superconducting state.
    type(superconductor) :: this

    ! Call the superclass constructor
    this%conductor = conductor()

    ! Initialize the order parameter
    allocate(this%gap(size(this%location)))
    call this%set_gap( (1.0_wp,0.0_wp) )

    ! Initialize the coupling constant
    if (this%energy(size(this%energy)) > 0) then
      this%coupling = 1/acosh(this%energy(size(this%energy)))
    end if
  end function

  pure subroutine superconductor_init(this, gap)
    ! Redefine the default initializer.
    class(superconductor), intent(inout) :: this
    complex(wp),           intent(in   ) :: gap

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
  end subroutine

  impure subroutine superconductor_update_gap(this)
    ! Update the superconducting gap if the BCS coupling constant is nonzero.
    class(superconductor), intent(inout) :: this         ! Superconductor object that will be updated
    complex(wp),           allocatable   :: gap(:)       ! Used to calculate the superconducting order parameter
    complex(wp)                          :: diff         ! Change in the superconducting order parameter mean
    complex(wp)                          :: singlet      ! Singlet component of the anomalous propagators
    integer                              :: n, m         ! Loop variables

    if (abs(this%coupling) > eps) then
      ! Allocate workspace memory
      allocate(gap(size(this%energy)))

      ! Calculate the mean superconducting order parameter
      diff = this%get_gap_mean()

      ! Iterate over the stored propagators
      do n = 1,size(this%location)
        do m = 1,size(this%energy)
          ! Calculate the singlet components of the anomalous propagators
          singlet = ( this%propagator(m,n)%get_f_s() - conjg(this%propagator(m,n)%get_ft_s()) )/2.0_wp

          ! Calculate the gap equation integrand and store it in an array
          gap(m)  = singlet * this%coupling * tanh(0.8819384944310228_wp * this%energy(m)/this%temperature)
        end do

        ! Interpolate and integrate the results, and update the superconducting order parameter
        this%gap(n) = integrate(this%energy, gap, 1e-6_wp, this%energy(ubound(this%energy,1)))
      end do

      ! Calculate the difference in mean superconducting order parameter
      diff = this%get_gap_mean() - diff

      ! Status information
      if (this%information >= 0 .and. .not. this % lock) then
        write(stdout,'(4x,a,f10.8,a)') 'Gap change: ',abs(diff),'                                        '
        flush(stdout)
      end if

      ! Deallocate workspace memory
      deallocate(gap)
    end if
  end subroutine

  !--------------------------------------------------------------------------------!
  !                    IMPLEMENTATION OF GETTERS AND SETTERS                       !
  !--------------------------------------------------------------------------------!

  pure subroutine superconductor_set_gap(this, gap)
    ! Updates the superconducting order parameter from a scalar.
    class(superconductor), intent(inout) :: this
    complex(wp),           intent(in   ) :: gap
    integer                              :: n

    do n = 1,size(this%gap)
      this%gap(n) = gap
    end do
  end subroutine

  pure function superconductor_get_gap(this, location) result(gap)
    ! Returns the superconducting order parameter at the given location.
    class(superconductor), intent(in) :: this
    real(wp),              intent(in) :: location
    complex(wp)                       :: gap
    integer                           :: n

    ! Calculate the index corresponding to the given location
    n = floor(location*(size(this%location)-1) + 1)

    ! Extract the superconducting order parameter at that point
    if (n <= 1) then
      gap = this%gap(1)
    else
      gap = this%gap(n-1) + (this%gap(n)-this%gap(n-1))*(location-this%location(n-1))/(this%location(n)-this%location(n-1))
    end if
  end function

  pure function superconductor_get_gap_mean(this) result(gap)
    ! Returns the superconducting order parameter average in the material.
    class(superconductor), intent(in)  :: this
    complex(wp)                        :: gap

    gap = sum(this%gap)/max(1,size(this%gap)) 
  end function

  !--------------------------------------------------------------------------------!
  !                      IMPLEMENTATION OF UTILITY METHODS                         !
  !--------------------------------------------------------------------------------!

  impure subroutine superconductor_conf(this, key, val)
    !! Configure a material property based on a key-value pair.
    class(superconductor), intent(inout) :: this
    character(*),          intent(in   ) :: key
    character(*),          intent(in   ) :: val

    select case(key)
      case("coupling")
        call evaluate(val, this%coupling)
      case ('gap')
        ! @TODO: Split this into two parts, namely 'gap' and 'phase'.
        !        * Constructor should initialize to (1,0).
        !        * Gap as a function of z, should change the magnitude:
        !            gap(n) -> gap(n)*(gap/abs(gap(n)))
        !        * Phase as a function of z, should change the phase:
        !            gap(n) -> abs(gap(n)) * exp((0,pi)*phase(n))
        !        * So gap and phase should become two distinct parameters.
        block
          real(wp) :: gap
          real(wp) :: phase
          integer  :: iostat

          iostat = 0
          read(val,*,iostat=iostat) gap, phase
          if ( iostat /= 0 ) then
            read(val,*) gap
            phase = 0
          end if

          call this % init( gap = gap * exp( (0.0,1.0)*pi*phase ) )
        end block
      case default
        call this%conductor%conf(key, val)
    end select
  end subroutine
end module
