!> Author:   Jabir Ali Ouassou
!> Date:     2015-07-20
!> Category: Materials
!>
!> This module defines the data type 'ferromagnet', which models the physical state of a ferromagnet. The type is a
!> member of class(conductor), and thus inherits the internal structure and generic methods defined in conductor_m.

module ferromagnet_m
  use stdio_m
  use math_m
  use spin_m
  use conductor_m
  implicit none
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

    ! These methods are used to access and mutate the parameters
    procedure                        :: set_exchange       => ferromagnet_set_exchange        ! Updates the ferromagnetic order parameter from a given vector

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
      if (norm2(this%exchange) < sqrt(eps)) then
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
          h(n)  = ((0.0_wp,-1.0_wp)/this%thouless) * (exchange(1,n)*pauli1 + exchange(2,n)*pauli2 + exchange(3,n)*pauli3)
          ht(n) = ((0.0_wp,+1.0_wp)/this%thouless) * (exchange(1,n)*pauli1 - exchange(2,n)*pauli2 + exchange(3,n)*pauli3)
        end do

        end associate
      end if
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
  !                    IMPLEMENTATION OF GETTERS AND SETTERS                       !
  !--------------------------------------------------------------------------------!

  pure subroutine ferromagnet_set_exchange(this, exchange)
    ! Updates the ferromagnetic order parameter from a vector.
    class(ferromagnet), intent(inout) :: this
    real(wp),           intent(in)    :: exchange(3)
    integer                           :: n

    ! Allocate the exchange field array if necessary
    if (.not. allocated(this%exchange)) then
      allocate(this%exchange(3,size(this%location)))
    end if

    ! Copy data from input array
    do n = 1,size(this%location)
      this%exchange(:,n) = exchange
    end do
  end subroutine

  !--------------------------------------------------------------------------------!
  !                      IMPLEMENTATION OF UTILITY METHODS                         !
  !--------------------------------------------------------------------------------!

  impure subroutine ferromagnet_conf(this, key, val)
    !! Configure a material property based on a key-value pair.
    class(ferromagnet), intent(inout) :: this
    character(*),       intent(in)    :: key
    character(*),       intent(in)    :: val

    select case(key)
      case ('magnetization')
        call evaluate(val, this % location, this % exchange)
      case default
        call this % conductor % conf(key, val)
    end select
  end subroutine
end module
