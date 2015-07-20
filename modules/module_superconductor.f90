! This module defines the data structure 'superconductor', which models the physical state of such
! materials as a function of position and energy. This data type inherits the internal structure of
! the 'conductor' type that was defined in module_conductor, and thus belongs to class(conductor).
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-17
! Updated: 2015-07-20

module module_superconductor
  use module_precision
  use module_spin
  use module_state
  use module_conductor
  implicit none

  ! Type declaration
  type, extends(conductor) :: superconductor
    real(dp)                 :: temperature = 1e-6_dp         ! Temperature of the system (relative to the critical temperature of a bulk superconductor)
    real(dp)                 :: coupling    = 0.20_dp               ! BCS coupling constant that defines the strength of the superconductor (dimensionless)
    complex(dp), allocatable :: gap(:)                             ! Superconducting gap as a function of position (relative to the zero-temperature gap of a bulk superconductor)
    contains
    procedure                :: get_gap          => superconductor_get_gap
    procedure                :: usadel_equation  => superconductor_usadel_equation
    procedure                :: update_fields    => superconductor_update_fields   
  end type

  ! Type constructor
  interface superconductor
    module procedure superconductor_construct_bcs
  end interface

contains
  pure function superconductor_construct_bcs(energy, gap, scattering, points) result(this)
    ! Constructs a superconductor object corresponding to a BCS superconductor with a given position and energy range
    type(superconductor)              :: this        ! Superconductor object that will be constructed
    real(dp),    intent(in)           :: energy(:)   ! Discretized energy domain that will be used
    real(dp),    intent(in), optional :: scattering  ! Imaginary energy term
    complex(dp), intent(in), optional :: gap         ! Superconducting gap 
    integer,     intent(in), optional :: points      ! Number of positions 

    ! Call the superclass constructor
    this%conductor = conductor_construct_bcs(energy, gap=gap, scattering=scattering, points=points)

    ! Allocate memory (if necessary)
    if (.not. allocated(this%gap)) then
      allocate(this%gap(size(this%conductor%location)))
    end if

    ! Initialize the superconducting gap
    if (present(gap)) then
      this%gap = gap
    else
      this%gap = (1.0_dp,0.0_dp)
    end if
  end function

  pure subroutine superconductor_destruct(this)
    ! Define the type destructor
    type(superconductor), intent(inout) :: this

    ! Deallocate memory (if necessary)
    if(allocated(this%gap)) then
      deallocate(this%gap)
    end if

    ! Call the superclass destructor
    call conductor_destruct(this%conductor)
  end subroutine

  pure function superconductor_get_gap(this, location) result(gap)
    ! Returns the superconducting order parameter at the given location.
    class(superconductor), intent(in) :: this
    real(dp),              intent(in) :: location
    complex(dp)                       :: gap
    integer                           :: n

    ! Calculate the index corresponding to the given location
    n = nint(location*(size(this%location)-1) + 1)

    ! Extract the superconducting order parameter at that point
    gap = this%gap(n)
  end function

  subroutine superconductor_usadel_equation(this, z, g, gt, dg, dgt, d2g, d2gt)
    ! Use the Usadel equation to calculate the second derivatives of the Riccati parameters at point z.
    class(superconductor), intent(in)  :: this
    real(dp),              intent(in)  :: z
    type(spin),            intent(in)  :: g, gt, dg, dgt
    type(spin),            intent(out) :: d2g, d2gt
    type(spin)                         :: N, Nt
    complex(dp)                        :: gap, gapt

    ! Lookup the superconducting order parameter
    gap  = this%get_gap(z)/this%thouless
    gapt = conjg(gap)

    ! Calculate the normalization matrices
    N   = spin_inv( pauli0 - g*gt )
    Nt  = spin_inv( pauli0 - gt*g )

    ! Calculate the second derivatives of the Riccati parameters
    d2g  = (-2.0_dp,0.0_dp)*dg*Nt*gt*dg - (0.0_dp,2.0_dp)*this%erg*g  - gap  * pauli2 + gapt * g*pauli2*g
    d2gt = (-2.0_dp,0.0_dp)*dgt*N*g*dgt - (0.0_dp,2.0_dp)*this%erg*gt + gapt * pauli2 - gap  * gt*pauli2*gt
  end subroutine

  subroutine superconductor_update_fields(this)
    class(superconductor), intent(inout) :: this

    real(dp), allocatable                :: gap_real(:), dgap_real(:)
    real(dp), allocatable                :: gap_imag(:), dgap_imag(:)
    complex(dp)                          :: singlet
    integer                              :: n, m, err
    real(dp), external                   :: dpchqa

    allocate(gap_real(size(this%energy)))
    allocate(gap_imag(size(this%energy)))
    allocate(dgap_real(size(this%energy)))
    allocate(dgap_imag(size(this%energy)))

    do n = 1,size(this%location)
      ! Calculate the real and imaginary parts of the gap equation integrand
      do m = 1,size(this%energy)
        singlet  = ( this%state(m,n)%get_f_s() - conjg(this%state(m,n)%get_ft_s()) )/2.0_dp

        gap_real(m) =  dble(singlet) * this%coupling * tanh(0.8819384944310228_dp * this%energy(m)/this%temperature)
        gap_imag(m) = aimag(singlet) * this%coupling * tanh(0.8819384944310228_dp * this%energy(m)/this%temperature)
      end do

      ! Create a Piecewise Cubic Hermitian Interpolation of the numerical results above
      call dpchez(size(this%energy), this%energy, gap_real, dgap_real, .false., 0, 0, err)
      call dpchez(size(this%energy), this%energy, gap_imag, dgap_imag, .false., 0, 0, err)

      ! Perform a numerical integration of the interpolation, and update the superconducting order parameter
      this%gap(n) = cmplx( dpchqa(size(this%energy), this%energy, gap_real, dgap_real, 0.0_dp, cosh(1.0_dp/this%coupling), err), &
                           dpchqa(size(this%energy), this%energy, gap_imag, dgap_imag, 0.0_dp, cosh(1.0_dp/this%coupling), err), &
                           kind=dp )
    end do

    deallocate(gap_real)
    deallocate(gap_imag)
    deallocate(dgap_real)
    deallocate(dgap_imag)
  end subroutine
end module
