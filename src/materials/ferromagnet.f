! This module defines the data type 'ferromagnet', which models the physical state of a ferromagnet. The type is a
! member of class(conductor), and thus inherits the internal structure and generic methods defined in mod_conductor.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-20
! Updated: 2015-08-10

module mod_ferromagnet
  use mod_conductor
  implicit none

  ! Type declaration
  type, extends(conductor)           :: ferromagnet
    real(dp),   allocatable          :: exchange(:,:)                                         ! Magnetic exchange field as a function of position
    type(spin), allocatable, private :: h(:), ht(:)                                           ! Used by internal subroutines to handle exchange fields
  contains
    procedure                        :: diffusion_equation => ferromagnet_diffusion_equation  ! Differential equation that describes the ferromagnet
    procedure                        :: update_prehook     => ferromagnet_update_prehook      ! Code to execute before calculating the Green's functions
    procedure                        :: update_posthook    => ferromagnet_update_posthook     ! Code to execute after  calculating the Green's functions
    final                            ::                       ferromagnet_destruct            ! Type destructor
  end type

  ! Type constructor
  interface ferromagnet
    module procedure ferromagnet_construct_homogeneous
  end interface
contains

  !--------------------------------------------------------------------------------!
  !                IMPLEMENTATION OF CONSTRUCTORS AND DESTRUCTORS                  !
  !--------------------------------------------------------------------------------!

  pure function ferromagnet_construct_homogeneous(energy, exchange, gap, thouless, scattering, points) result(this)
    ! Constructs a ferromagnet object initialized to a weak superconductor.
    type(ferromagnet)                 :: this         ! Ferromagnet object that will be constructed
    real(dp),    intent(in)           :: energy(:)    ! Discretized energy domain that will be used
    real(dp),    intent(in), optional :: exchange(3)  ! Magnetic exchange field
    complex(dp), intent(in), optional :: gap          ! Superconducting gap   (default: conductor default)
    real(dp),    intent(in), optional :: thouless     ! Thouless energy       (default: conductor default)
    real(dp),    intent(in), optional :: scattering   ! Imaginary energy term (default: conductor default)
    integer,     intent(in), optional :: points       ! Number of positions   (default: conductor default)
    integer                           :: n            ! Loop variable

    ! Call the superclass constructor
    this%conductor = conductor_construct(energy, gap=gap, thouless=thouless, scattering=scattering, points=points)


    ! Handle the exchange field argument
    if (present(exchange) .and. norm2(exchange) > 1e-10) then
      ! Allocate the exchange field array if necessary
      if (.not. allocated(this%exchange)) then
        allocate(this%exchange(3,size(this%location)))
      end if

      ! Copy data from input array
      do n = 1,size(this%location)
        this%exchange(:,n) = exchange
      end do
    end if
  end function

  pure subroutine ferromagnet_destruct(this)
    ! Define the type destructor.
    type(ferromagnet), intent(inout) :: this

    ! Deallocate memory (if necessary)
    if (allocated(this%exchange)) then
      deallocate(this%exchange)
    end if

    if (allocated(this%h)) then
      deallocate(this%h)
      deallocate(this%ht)
    end if

    ! Call the superclass destructor
    call conductor_destruct(this%conductor)
  end subroutine

  !--------------------------------------------------------------------------------!
  !                    IMPLEMENTATION OF FERROMAGNET METHODS                       !
  !--------------------------------------------------------------------------------!

  pure subroutine ferromagnet_diffusion_equation(this, e, z, g, gt, dg, dgt, d2g, d2gt)
    ! Use the diffusion equation to calculate the second derivatives of the Riccati parameters at point z.
    class(ferromagnet), intent(in   ) :: this
    complex(dp),        intent(in   ) :: e
    real(dp),           intent(in   ) :: z
    type(spin),         intent(in   ) :: g, gt, dg, dgt
    type(spin),         intent(inout) :: d2g, d2gt
    type(spin)                        :: h, ht
    real(dp)                          :: d
    integer                           :: n

    ! Calculate the second derivatives of the Riccati parameters (conductor terms)
    call this%conductor%diffusion_equation(e, z, g, gt, dg, dgt, d2g, d2gt)

    if (allocated(this%exchange)) then
      ! Calculate the index corresponding to the given location
      n = floor(z*(size(this%location)-1) + 1)

      ! Extract the exchange field terms at that point
      if (n <= 1) then
        h  = this%h(1)
        ht = this%ht(1)
      else
        d  = (z - this%location(n-1))/(this%location(n) - this%location(n-1))
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
      if (norm2(this%exchange) < 1e-10) then
        ! Deallocate negligible fields
        deallocate(this%exchange)
      else
        ! Allocate space for internal variables
        if (.not. allocated(this%h)) then
          allocate(this%h(size(this%location)))
          allocate(this%ht(size(this%location)))
        end if

        ! Rename the internal variables
        associate(exchange => this % exchange, &
                  h        => this % h,        &
                  ht       => this % ht)

        ! Update the internal variables
        do n = 1,size(exchange,2)
          h(n)  = ((0.0_dp,-1.0_dp)/this%thouless) * (exchange(1,n)*pauli1 + exchange(2,n)*pauli2 + exchange(3,n)*pauli3)
          ht(n) = ((0.0_dp,+1.0_dp)/this%thouless) * (exchange(1,n)*pauli1 - exchange(2,n)*pauli2 + exchange(3,n)*pauli3)
        end do

        end associate
      end if
    end if

    ! Modify the type string
    if (allocated(this%exchange)) then
      if (colors) then
        this%type_string = color_red // 'FERROMAGNET' // color_none
        if (allocated(this%spinorbit))       this%type_string = trim(this%type_string) // color_cyan   // ' [SOC]' // color_none
        if (allocated(this%magnetization_a)) this%type_string = trim(this%type_string) // color_purple // ' [SAL]' // color_none
        if (allocated(this%magnetization_b)) this%type_string = trim(this%type_string) // color_purple // ' [SAR]' // color_none
      else
        this%type_string = 'FERROMAGNET'
        this%type_string = 'SUPERCONDUCTOR'
        if (allocated(this%spinorbit))       this%type_string = trim(this%type_string) // ' [SOC]'
        if (allocated(this%magnetization_a)) this%type_string = trim(this%type_string) // ' [SAL]'
        if (allocated(this%magnetization_b)) this%type_string = trim(this%type_string) // ' [SAR]'
      end if
    end if
  end subroutine

  impure subroutine ferromagnet_update_posthook(this)
    ! Code to execute after running the update method of a class(ferromagnet) object.
    class(ferromagnet), intent(inout) :: this

    continue
  end subroutine
end module
