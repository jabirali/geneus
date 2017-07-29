!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This submodule is included by conductor.f, and contains the equations which model spin-flip and spin-orbit scattering.
!>
!> @TODO 
!>   Rewrite using the new nambu.f library, replacing e.g. diag(m·σ,m·σ*) with m*nambuv(1:3).
!>   Also, we may then for brevity replace matmul(G,matmul(M,G)) with G*M*G, and so on.
!>
!> @TODO
!>   Move the orbital depairing implementation here. It might run slightly slower that way,
!>   but will be easier to maintain, especially when generalizing to nonequilibrium.


module spinscattering_m
  use :: material_m
  use :: math_m
  use :: spin_m
  use :: propagator_m
  private

  ! Public interface
  public spinscattering, spinscattering_construct

  ! Type declarations
  type :: spinscattering
    class(material), pointer    :: material  => null()                   !! Pointer to the material modelled by this instance
    real(wp)                    :: spinflip  =  0.0_wp                   !! Spin-flip  scattering coefficient (definition: 1/8Δτ)
    real(wp)                    :: spinorbit =  0.0_wp                   !! Spin-orbit scattering coefficient (definition: 1/8Δτ)
    complex(wp), dimension(4,4) :: sigma1, sigma2, sigma3, tau3          !! Pauli matrices that are used internally in this object
  contains
    procedure :: diffusion_equation => spinscattering_diffusion_equation !! Defines the Usadel diffusion equation (scattering terms)
  end type

  ! Type constructors
  interface spinscattering
    module procedure spinscattering_construct
  end interface
contains
  pure subroutine spinscattering_diffusion_equation(this, g, gt, dg, dgt, d2g, d2gt) 
    !! Calculate the spin-flip and spin-orbit scattering terms in the diffusion
    !! equation, and update the second derivatives of the Riccati parameters.
    class(spinscattering), intent(in)         :: this
    type(spin),            intent(in)         :: g, gt, dg, dgt
    type(spin),            intent(inout)      :: d2g, d2gt
    complex(wp),           dimension(4,4)     :: p, u, v
    real(wp)                                  :: sf, so
    type(propagator)                          :: GR

    ! Calculate the 4×4 propagator matrix
    GR = propagator(g, gt)
    p  = GR % retarded()

    ! Calculate the spin-flip and spin-orbit coefficients
    sf = this%spinflip  / (2 * this%material%thouless)
    so = this%spinorbit / (2 * this%material%thouless)

    ! Calculate the 4×4 diffusion equation contribution
    associate(s1 => this % sigma1, s2 => this % sigma2, s3 => this % sigma3, t3 => this % tau3)
      v = matmul(s1,matmul(p,s1)) + matmul(s2,matmul(p,s2)) + matmul(s3,matmul(p,s3))
      u = commutator(sf*v + so*matmul(t3,matmul(v,t3)), p)
    end associate

    ! Update the second derivatives of the Riccati parameters
    d2g  = d2g  + (pauli0 - g*gt) * (u(1:2,3:4) - u(1:2,1:2)*g )
    d2gt = d2gt + (pauli0 - gt*g) * (u(3:4,1:2) - u(3:4,3:4)*gt)
  end subroutine

  function spinscattering_construct(parent) result(this)
    !! Constructs a spinscattering object with a given parent material.
    type(spinscattering)    :: this
    class(material), target :: parent

    ! Save a pointer to the parent object
    this % material => parent

    ! Define the necessary 4×4 basis matrices
    this % tau3   = diag(+pauli0%matrix, -pauli0%matrix)
    this % sigma1 = diag(+pauli1%matrix, +pauli1%matrix)
    this % sigma2 = diag(+pauli2%matrix, -pauli2%matrix)
    this % sigma3 = diag(+pauli3%matrix, +pauli3%matrix)
  end function
end module
