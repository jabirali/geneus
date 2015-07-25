! This module defines the data type 'ferromagnet', which models the physical state of a ferromagnet. The type is a
! member of class(conductor), and thus inherits the internal structure and generic methods defined in mod_conductor.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-20
! Updated: 2015-07-25

module mod_ferromagnet
  use mod_system
  use mod_spin
  use mod_green
  use mod_conductor
  implicit none

  ! Type declaration
  type, extends(conductor) :: ferromagnet
    real(dp), allocatable    :: exchange(:,:)                                         ! Magnetic exchange field as a function of position
  contains
    procedure                :: usadel_equation       => ferromagnet_usadel_equation  ! Differential equation that describes the superconductor
    procedure                :: get_exchange          => ferromagnet_get_exchange     ! Returns the magnetic exchange field at a given position
    final                    ::                          ferromagnet_destruct         ! Type destructor
  end type

  ! Type constructor
  interface ferromagnet
    module procedure ferromagnet_construct_homogeneous
  end interface
contains

  !--------------------------------------------------------------------------------!
  !                IMPLEMENTATION OF CONSTRUCTORS AND DESTRUCTORS                  !
  !--------------------------------------------------------------------------------!

  pure function ferromagnet_construct_homogeneous(energy, exchange, gap, thouless, scattering, points, spinorbit) result(this)
    ! Constructs a ferromagnet object initialized to a weak superconductor.
    type(ferromagnet)                 :: this         ! Ferromagnet object that will be constructed
    real(dp),    intent(in)           :: energy(:)    ! Discretized energy domain that will be used
    real(dp),    intent(in)           :: exchange(3)  ! Magnetic exchange field
    complex(dp), intent(in), optional :: gap          ! Superconducting gap   (default: conductor default)
    real(dp),    intent(in), optional :: thouless     ! Thouless energy       (default: conductor default)
    real(dp),    intent(in), optional :: scattering   ! Imaginary energy term (default: conductor default)
    integer,     intent(in), optional :: points       ! Number of positions   (default: conductor default)
    type(spin),  intent(in), optional :: spinorbit(:) ! Spin-orbit coupling   (default: zero coupling)
    integer                           :: n            ! Loop variable

    ! Call the superclass constructor
    this%conductor = conductor_construct(energy, gap=gap, thouless=thouless, scattering=scattering, &
                                         points=points, spinorbit=spinorbit)

    ! Allocate memory (if necessary)
    if (.not. allocated(this%exchange)) then
      allocate(this%exchange(3,size(this%conductor%location)))
    end if

    ! Initialize the exchange field
    do n = 1,size(this%conductor%location)
      this%exchange(:,n) = exchange
    end do

    ! Modify the type string
    this%type_string = color_red // 'FERROMAGNET' // color_none
    if (allocated(this%spinorbit)) then
      this%type_string = trim(this%type_string) // color_cyan // ' [SOC] ' // color_none
    end if
  end function

  pure subroutine ferromagnet_destruct(this)
    ! Define the type destructor.
    type(ferromagnet), intent(inout) :: this

    ! Deallocate memory (if necessary)
    if(allocated(this%exchange)) then
      deallocate(this%exchange)
    end if

    ! Call the superclass destructor
    call conductor_destruct(this%conductor)
  end subroutine

  !--------------------------------------------------------------------------------!
  !                    IMPLEMENTATION OF FERROMAGNET METHODS                       !
  !--------------------------------------------------------------------------------!

  subroutine ferromagnet_usadel_equation(this, z, g, gt, dg, dgt, d2g, d2gt)
    ! Use the Usadel equation to calculate the second derivatives of the Riccati parameters at point z.
    class(ferromagnet), intent(in   ) :: this
    real(dp),           intent(in   ) :: z
    type(spin),         intent(in   ) :: g, gt, dg, dgt
    type(spin),         intent(inout) :: d2g, d2gt
    type(spin)                        :: P, Pt
    real(dp)                          :: h(3)

    ! Lookup the magnetic exchange field
    h  = this%get_exchange(z)/this%thouless

    ! Calculate the exchange field factors in the Usadel equation
    P  = (0.0_dp,-1.0_dp) * (h(1)*pauli1 + h(2)*pauli2 + h(3)*pauli3)
    Pt = (0.0_dp,+1.0_dp) * (h(1)*pauli1 - h(2)*pauli2 + h(3)*pauli3)

    ! Calculate the second derivatives of the Riccati parameters (conductor terms)
    call this%conductor%usadel_equation(z, g, gt, dg, dgt, d2g, d2gt)

    ! Calculate the second derivatives of the Riccati parameters (ferromagnet terms)
    d2g  = d2g  + P  * g  + g  * Pt
    d2gt = d2gt + Pt * gt + gt * P
  end subroutine

  !--------------------------------------------------------------------------------!
  !                    IMPLEMENTATION OF GETTERS AND SETTERS                       !
  !--------------------------------------------------------------------------------!

  pure function ferromagnet_get_exchange(this, location) result(h)
    ! Returns the magnetic exchange field at the given location.
    class(ferromagnet), intent(in) :: this
    real(dp),           intent(in) :: location
    real(dp)                       :: h(3)
    integer                        :: n

    ! Calculate the index corresponding to the given location
    n = nint(location*(size(this%location)-1) + 1)

    ! Extract the magnetic exchange field at that point
    h = this%exchange(:,n)
  end function
end module
