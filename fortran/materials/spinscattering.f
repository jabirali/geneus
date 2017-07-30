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
    type(nambu)                          :: G, R
    type(nambu)                          :: Gsf, Gso, Gdp
    real(wp)                             :: Csf, Cso, Cdp

    ! Construct the propagator matrix
    G = p % retarded()

    ! Construct the self-energy matrices
    associate(N => this % nambuv)
      Gdp = N(4)*G*N(4)
      Gsf = N(1)*G*N(1) + N(2)*G*N(2) + N(3)*G*N(3)
      Gso = N(4)*Gsf*N(4)
    end associate

    ! Calculate the self-energy prefactors
    Cdp = (this % depairing) / (8 * this % material % thouless)
    Csf = (this % spinflip ) / (2 * this % material % thouless)
    Cso = (this % spinorbit) / (2 * this % material % thouless)

    ! Calculate the self-energy commutators
    R = Cdp*Gdp + Csf*Gsf + Cso*Gso
    R = R*G - G*R

    ! Update the second derivatives of the Riccati parameters
    associate(  U => R % matrix,                  &
                g => p % g,       gt => p % gt,   &
               dg => p % dg,     dgt => p % dgt,  &
              d2g => p % d2g,   d2gt => p % d2gt  )
      d2g  = d2g  + (pauli0 - g*gt) * (U(1:2,3:4) - U(1:2,1:2)*g )
      d2gt = d2gt + (pauli0 - gt*g) * (U(3:4,1:2) - U(3:4,3:4)*gt)
    end associate
  end subroutine
end module
