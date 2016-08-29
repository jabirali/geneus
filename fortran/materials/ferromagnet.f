!> Author:   Jabir Ali Ouassou
!> Date:     2015-07-20
!> Category: Materials
!>
!> This module defines the data type 'ferromagnet', which models the physical state of a ferromagnet. The type is a
!> member of class(conductor), and thus inherits the internal structure and generic methods defined in conductor_m.

module ferromagnet_m
  use :: stdio_m
  use :: math_m
  use :: spin_m
  use :: conductor_m
  private

  ! Type declaration
  type, public, extends(conductor)   :: ferromagnet
    real(wp),   allocatable          :: exchange(:,:)                                         ! Magnetic exchange field as a function of position
    type(spin), allocatable, private :: h(:), ht(:)                                           ! Used by internal subroutines to handle exchange fields
  contains
    ! These methods contain the equations that describe ferromagnets
    procedure                        :: diffusion_equation => ferromagnet_diffusion_equation  ! Differential equation that describes the ferromagnet
    procedure                        :: update_prehook     => ferromagnet_update_prehook      ! Code to execute before calculating the propagators
    procedure                        :: update_posthook    => ferromagnet_update_posthook     ! Code to execute after  calculating the propagators

    ! These methods define miscellaneous utility functions
    procedure                        :: conf               => ferromagnet_conf                ! Configures material parameters
  end type

  ! Type constructor
  interface ferromagnet
    module procedure ferromagnet_construct
  end interface
contains

  !--------------------------------------------------------------------------------!
  !                        IMPLEMENTATION OF CONSTRUCTORS                          !
  !--------------------------------------------------------------------------------!

  function ferromagnet_construct() result(this)
    ! Constructs a ferromagnetic material that is initialized as a weak superconductor.
    type(ferromagnet) :: this

    ! Call the superclass constructor
    this%conductor = conductor()
  end function

  !--------------------------------------------------------------------------------!
  !                    IMPLEMENTATION OF FERROMAGNET METHODS                       !
  !--------------------------------------------------------------------------------!

  pure subroutine ferromagnet_diffusion_equation(this, e, z, g, gt, dg, dgt, d2g, d2gt)
    ! Use the diffusion equation to calculate the second derivatives of the Riccati parameters at point z.
    class(ferromagnet), intent(in)    :: this
    complex(wp),        intent(in)    :: e
    real(wp),           intent(in)    :: z
    type(spin),         intent(in)    :: g, gt, dg, dgt
    type(spin),         intent(inout) :: d2g, d2gt
    type(spin)                        :: h, ht
    real(wp)                          :: d
    integer                           :: n, m

    ! Calculate the second derivatives of the Riccati parameters (conductor terms)
    call this%conductor%diffusion_equation(e, z, g, gt, dg, dgt, d2g, d2gt)

    if (allocated(this%h)) then
      ! Calculate the index corresponding to the given location
      m = size(this%h) - 1          ! Number of array intervals
      n = floor(z*m + 1)            ! Nearest position in array

      ! Extract the exchange field terms at that point
      if (n <= 1) then
        ! Left edge of the material
        h  = this%h(1)
        ht = this%ht(1)
      else
        ! Linear interpolation from known values. The relative displacement d is defined
        ! as [z - location(n-1)]/[location(n) - location(n-1)], but assuming location(:)
        ! is a uniform array of values in the range [0,1], the below will be equivalent.
        d  = z*m - (n-2)
        h  = this%h(n-1)  + (this%h(n)  - this%h(n-1))  * d
        ht = this%ht(n-1) + (this%ht(n) - this%ht(n-1)) * d
      end if

      ! Calculate the second derivatives of the Riccati parameters (ferromagnet terms)
      d2g  = d2g  + h  * g  + g  * ht
      d2gt = d2gt + ht * gt + gt * h
    end if
  end subroutine

  impure subroutine ferromagnet_update_prehook(this)
    ! Updates the exchange field terms in the diffusion equation.
    class(ferromagnet), intent(inout) :: this ! Ferromagnet object that will be updated
    integer                           :: n    ! Loop variable

    ! Call the superclass prehook
    call this%conductor%update_prehook

    ! Rename the internal variables
    if (allocated(this%exchange)) then
      ! Allocate space for internal variables
      if (.not. allocated(this%h)) then
        allocate(this%h (size(this%exchange,2)))
        allocate(this%ht(size(this%exchange,2)))
      end if

      ! Update the internal variables
      associate(exchange => this % exchange, h => this % h, ht => this % ht)
        do n = 1,size(exchange,2)
          h(n)  = ((0.0_wp,-1.0_wp)/this%thouless) * (exchange(1,n)*pauli1 + exchange(2,n)*pauli2 + exchange(3,n)*pauli3)
          ht(n) = ((0.0_wp,+1.0_wp)/this%thouless) * (exchange(1,n)*pauli1 - exchange(2,n)*pauli2 + exchange(3,n)*pauli3)
        end do
      end associate
    end if

    ! Modify the type string
    if (allocated(this%exchange)) then
      this%type_string = color_red // 'FERROMAGNET' // color_none
      if (allocated(this%spinorbit))       this%type_string = trim(this%type_string) // color_cyan   // ' [SOC]' // color_none
      if (allocated(this%magnetization_a)) this%type_string = trim(this%type_string) // color_purple // ' [SAL]' // color_none
      if (allocated(this%magnetization_b)) this%type_string = trim(this%type_string) // color_purple // ' [SAR]' // color_none
    end if
  end subroutine

  impure subroutine ferromagnet_update_posthook(this)
    ! Code to execute after running the update method of a class(ferromagnet) object.
    class(ferromagnet), intent(inout) :: this

    ! Call the superclass posthook
    call this%conductor%update_posthook
  end subroutine

  !--------------------------------------------------------------------------------!
  !                      IMPLEMENTATION OF UTILITY METHODS                         !
  !--------------------------------------------------------------------------------!

  impure subroutine ferromagnet_conf(this, key, val)
    !! Configure a material property based on a key-value pair.
    use :: calculus_m
    use :: evaluate_m

    class(ferromagnet), intent(inout)   :: this
    character(*),       intent(in)      :: key
    character(*),       intent(in)      :: val

    select case(key)
      case ('magnetization')
        block
          ! Allocate memory for the location array
          real(wp), allocatable, dimension(:) :: location
          allocate(location(1024 * size(this%location)))

          ! Discretize the locations in the material
          call linspace(location, 0.0_wp, 1.0_wp)

          ! Initialize the magnetic exchange field
          call evaluate(val, location, this % exchange)

          ! Deallocate the field if it is negligible
          if (norm2(this % exchange) < sqrt(eps)) then
            deallocate(this % exchange)
          end if

          ! Deallocate the location array
          deallocate(location)
        end block
      case default
        ! Pass this option to the superclass
        call this % conductor % conf(key, val)
    end select
  end subroutine
end module
