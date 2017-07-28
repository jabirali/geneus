!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This submodule is included by conductor.f, and contains the equations which model spin-active interfaces.
!>
!> @TODO: Rewrite using the new nambu.f library, replacing e.g. diag(m·σ,m·σ*) with m*nambuv(1:3).
!>        Also, we may then for brevity replace matmul(G,matmul(M,G)) with G*M*G, and so on.

module spinactive_m
  use :: propagator_m
  use :: material_m
  use :: math_m
  use :: spin_m
  use :: nambu_m

  ! Type declarations
  type, public :: spinactive
    real(wp)                 :: conductance   = 0.0      !! Interfacial conductance
    real(wp)                 :: polarization  = 0.0      !! Interfacial spin-polarization
    real(wp)                 :: spinmixing    = 0.0      !! Interfacial 1st-order spin-mixing
    real(wp)                 :: secondorder   = 0.0      !! Interfacial 2nd-order spin-mixing
    real(wp), dimension(1:3) :: magnetization = [0,0,1]  !! Interfacial magnetization direction
    real(wp), dimension(1:3) :: misalignment  = [0,0,0]  !! Interfacial magnetization misalignment
  
    type(nambu), private     :: M                        !! Magnetization matrix (transmission)
    type(nambu), private     :: M0                       !! Magnetization matrix (reflection, this  side)
    type(nambu), private     :: M1                       !! Magnetization matrix (reflection, other side)
  contains
    procedure :: diffusion_equation_a => spinactive_diffusion_equation_a
    procedure :: diffusion_equation_b => spinactive_diffusion_equation_b
    procedure :: update_prehook       => spinactive_update_prehook
  end type
contains
  pure subroutine spinactive_update_prehook(this)
    !! Updates the internal variables associated with spin-active interfaces.
    class(spinactive), intent(inout) :: this 
  
    ! Process transmission properties
    this % M = nambuv(this % magnetization)
  
    ! Default reflection properties match transmission properties
    this % M0 = this % M
    this % M1 = this % M
  
    ! Process reflection properties (this side)
    if (nonzero(this % misalignment)) then
      this % M0 = nambuv(this % misalignment)
    end if
  
    ! ! Process reflection properties (other side)
    ! if (associated(this % material_a)) then
    !   select type (other => this % material_a)
    !     class is (conductor)
    !       call update_magnetization(this % spinactive_a % M1, other % spinactive_b % misalignment)
    !   end select
    ! end if
  
    ! if (associated(this % material_b)) then
    !   select type (other => this % material_b)
    !     class is (conductor)
    !       call update_magnetization(this % spinactive_b % M1, other % spinactive_a % misalignment)
    !   end select
    ! end if
  end subroutine
  
  pure subroutine spinactive_diffusion_equation_a(this, a, g1, gt1, dg1, dgt1, r1, rt1)
    !! Calculate the spin-active terms in the left boundary condition, and update the residuals.
    class(spinactive), target,  intent(in)    :: this
    type(propagator),           intent(in)    :: a
    type(spin),                 intent(in)    :: g1, gt1, dg1, dgt1
    type(spin),                 intent(inout) :: r1, rt1
    type(nambu)                               :: GM0, GM1, I
    type(propagator)                          :: GR0, GR1
  
    ! Calculate the 4×4 matrix propagators
    associate(g0 => a % g, gt0 => a % gt)
      GR0 = propagator(g0, gt0)
      GR1 = propagator(g1, gt1)

      GM0 = GR0 % retarded()
      GM1 = GR1 % retarded()
    end associate
  
    ! Calculate the 4×4 matrix current
    I = 0.25 * this%conductance &
      * spinactive_current(GM1, GM0, this%M, this%M0, this%M1, &
        this%polarization, this%spinmixing, this%secondorder)
  
    ! Calculate the deviation from the boundary condition
    r1  = dg1  + (pauli0 - g1*gt1) * (I % matrix(1:2,3:4) - I % matrix(1:2,1:2)*g1)
    rt1 = dgt1 + (pauli0 - gt1*g1) * (I % matrix(3:4,1:2) - I % matrix(3:4,3:4)*gt1)
  end subroutine
  
  pure subroutine spinactive_diffusion_equation_b(this, b, g2, gt2, dg2, dgt2, r2, rt2)
    !! Calculate the spin-active terms in the right boundary condition, and update the residuals.
    class(spinactive), target,  intent(in)    :: this
    type(propagator),           intent(in)    :: b
    type(spin),                 intent(in)    :: g2, gt2, dg2, dgt2
    type(spin),                 intent(inout) :: r2, rt2
    type(nambu)                               :: GM2, GM3, I
    type(propagator)                          :: GR2, GR3
  
    ! Calculate the 4×4 matrix propagators
    associate(g3 => b % g, gt3 => b % gt)
      GR2 = propagator(g2, gt2)
      GR3 = propagator(g3, gt3)

      GM2 = GR2 % retarded()
      GM3 = GR3 % retarded()
    end associate
  
    ! Calculate the 4×4 matrix current
    I = 0.25 * this%conductance &
      * spinactive_current(GM2, GM3, this%M, this%M0, this%M1, &
      this%polarization, this%spinmixing, this%secondorder)
  
    ! Calculate the deviation from the boundary condition
    r2  = dg2  - (pauli0 - g2*gt2) * (I % matrix(1:2,3:4) - I % matrix(1:2,1:2)*g2)
    rt2 = dgt2 - (pauli0 - gt2*g2) * (I % matrix(3:4,1:2) - I % matrix(3:4,3:4)*gt2)
  end subroutine
  
  pure function spinactive_current(G0, G1, M, M0, M1, P, Q, R) result(I)
    !! Calculate the matrix current at an interface with spin-active properties. The equations
    !! implemented here should be valid for an arbitrary interface polarization, and up to 2nd
    !! order in the transmission probabilities and spin-mixing angles of the interface. 
    type(nambu), intent(in) :: G0, G1      !! Propagator matrices
    type(nambu), intent(in) :: M0, M1, M   !! Magnetization matrices 
    real(wp),    intent(in) :: P,  Q,  R   !! Interface parameters
    type(nambu)             :: S0, S1      !! Matrix expressions
    type(nambu)             :: I           !! Matrix current
  
    ! Shortcut-evaluation for nonmagnetic interfaces
    if (abs(Q) == 0 .and. abs(P) == 0) then
      I = (G0*G1 - G1*G0)
      return
    end if
  
    ! Evaluate the 1st-order matrix functions
    S0 = spinactive_current1_transmission(G1)
    S1 = spinactive_current1_reflection()
  
    ! Evaluate the 1st-order matrix current
    associate(S => S0 + S1)
      I  = (G0*S - S*G0)
    end associate
  
    ! Calculate the 2nd-order contributions to the matrix current. Note that we make a
    ! number of simplifications in this implementation. In particular, we assume that 
    ! all interface parameters except the magnetization directions are equal on both
    ! sides of the interface. We also assume that the spin-mixing angles and tunneling
    ! probabilities of different channels have standard deviations that are much smaller
    ! than their mean values, which reduces the number of new fitting parameters to one.
  
    if (abs(R) > 0) then
      ! Evaluate the 1st-order matrix functions
      S1 = spinactive_current1_transmission(G1*M1*G1 - M1) 
  
      ! Evaluate the 2nd-order matrix current
      I = I                                     &
        + spinactive_current2_transmission()    &
        + spinactive_current2_crossterms()      &
        + spinactive_current2_reflection()
    end if
  contains
    pure function spinactive_current1_transmission(G) result(F)
      !! Calculate the 1st-order transmission terms in the matrix current commutator.
      type(nambu), intent(in) :: G
      type(nambu)             :: F
      real(wp) :: Pr, Pp, Pm
  
      Pr = sqrt(1 - P**2)
      Pp = 1 + Pr
      Pm = 1 - Pr
  
      F = G + (P/Pp)*(M*G+G*M) + (Pm/Pp)*(M*G*M)
    end function
  
    pure function spinactive_current1_reflection() result(F)
      !! Calculate the 1st-order spin-mixing terms in the matrix current commutator.
      type(nambu) :: F
  
      F = ((0,-1)*Q) * M0
    end function
  
    pure function spinactive_current2_transmission() result(I)
      !! Calculate the 2nd-order transmission terms in the matrix current.
      type(nambu) :: I
  
      I = (0.50*R/Q) * (S0*G0*S0)
    end function
  
    pure function spinactive_current2_reflection() result(I)
      !! Calculate the 2nd-order spin-mixing terms in the matrix current.
      type(nambu) :: I
  
      associate(U => M0*G0*M0)
        I = (0.25*R*Q) * (G0*U - U*G0)
      end associate
    end function
  
    pure function spinactive_current2_crossterms() result(I)
      !! Calculate the 2nd-order cross-terms in the matrix current.
      type(nambu) :: I
  
      associate(U => S0*G0*M0 + M0*G0*S0 + S1)
        I = ((0.00,0.25)*R) * (G0*U - U*G0)
      end associate
    end function
  end function
end module
