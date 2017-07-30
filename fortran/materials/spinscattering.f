!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This submodule is included by conductor.f, and contains equations which model 
!> spin-flip scattering, spin-orbit scattering, and magnetic orbital depairing.

module spinscattering_m
  use :: material_m
  use :: math_m
  use :: spin_m
  use :: nambu_m
  use :: propagator_m
  private

  ! Public interface
  public spinscattering, spinscattering_construct

  ! Type declarations
  type :: spinscattering
    ! Metadata
    class(material), pointer    :: material  => null()                     !! Pointer to the material modelled by this instance
    type(nambu), dimension(0:7) :: nambuv                                  !! Pauli matrices spanning the 4×4 Spin-Nambu space

    ! Physical fields
    real(wp)                    :: depairing =  0.0_wp                     !! Orbital depairing coefficient
    real(wp)                    :: spinflip  =  0.0_wp                     !! Spin-flip  scattering coefficient (1/8τΔ)
    real(wp)                    :: spinorbit =  0.0_wp                     !! Spin-orbit scattering coefficient (1/8τΔ)
  contains
    procedure :: diffusion_equation => spinscattering_diffusion_equation   !! Diffusion equation
  end type

  ! Type constructors
  interface spinscattering
    module procedure spinscattering_construct
  end interface
contains
  impure function spinscattering_construct(parent) result(this)
    !! Constructs a spinscattering object with a given parent material.
    type(spinscattering)    :: this
    class(material), target :: parent
    integer                 :: i

    ! Save a pointer to the parent object
    this % material => parent

    ! Memoize the most used basis matrices
    do i=0,7
      this % nambuv(i) = nambuv(i)
    end do
  end function

  pure subroutine spinscattering_diffusion_equation(this, p)
    !! Calculate the spin-flip and spin-orbit scattering terms in the diffusion
    !! equation, and update the second derivatives of the Riccati parameters.
    class(spinscattering), intent(in)    :: this
    type(propagator),      intent(inout) :: p
    type(nambu)                          :: g, r
    type(nambu)                          :: gsf, gso, gdp
    real(wp)                             :: csf, cso, cdp

    ! Construct the propagator matrix
    g = p % retarded()

    ! Construct the self-energy matrices
    associate(n => this % nambuv)
      gdp = n(4)*g*n(4)
      gsf = n(1)*g*n(1) + n(2)*g*n(2) + n(3)*g*n(3)
      gso = n(4)*gsf*n(4)
    end associate

    ! Calculate the self-energy prefactors
    cdp = (this % depairing) / (8 * this % material % thouless)
    csf = (this % spinflip ) / (2 * this % material % thouless)
    cso = (this % spinorbit) / (2 * this % material % thouless)

    ! Calculate the self-energy commutators
    r = cdp*gdp + csf*gsf + cso*gso
    r = r*g - g*r

    ! Update the second derivatives of the Riccati parameters
    associate(  u => r % matrix,                  &
                g => p % g,       gt => p % gt,   &
               dg => p % dg,     dgt => p % dgt,  &
              d2g => p % d2g,   d2gt => p % d2gt  )
      d2g  = d2g  + (pauli0 - g*gt) * (u(1:2,3:4) - u(1:2,1:2)*g )
      d2gt = d2gt + (pauli0 - gt*g) * (u(3:4,1:2) - u(3:4,3:4)*gt)
    end associate
  end subroutine
end module
