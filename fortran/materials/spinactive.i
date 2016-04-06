! This submodule is included by conductor.f, and contains the equations which model spin-active tunneling and vacuum interfaces.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-10-01
! Updated: 2016-04-06

! @TODO: Double-check what to do at vacuum interfaces with spinmixing.

pure subroutine spinactive_update_prehook(this)
  ! Updates the internal variables associated with spin-active interfaces.
  class(conductor), intent(inout) :: this 

  if (allocated(this % magnetization_a) .and. .not. this % reflecting_a) then
    ! Calculate the elements of the magnetization matrix diag(m·σ,m·σ*)
    associate(H => this % magnetization_a, M => this % M_a)
      M(1:2,1:2) = H(1) * pauli1 + H(2) * pauli2 + H(3) * pauli3
      M(3:4,3:4) = H(1) * pauli1 - H(2) * pauli2 + H(3) * pauli3
    end associate
  end if

  if (allocated(this % magnetization_b) .and. .not. this % reflecting_b) then
    ! Calculate the elements of the magnetization matrix diag(m·σ,m·σ*)
    associate(H => this % magnetization_b, M => this % M_b)
      M(1:2,1:2) = H(1) * pauli1 + H(2) * pauli2 + H(3) * pauli3
      M(3:4,3:4) = H(1) * pauli1 - H(2) * pauli2 + H(3) * pauli3
    end associate
  end if
end subroutine

pure subroutine spinactive_interface_equation_a(this, a, g1, gt1, dg1, dgt1, r1, rt1)
  ! Calculate the spin-active terms in the left boundary condition, and update the residuals.
  class(conductor), target, intent(in   ) :: this
  type(green),              intent(in   ) :: a
  type(spin),               intent(in   ) :: g1, gt1, dg1, dgt1
  type(spin),               intent(inout) :: r1, rt1
  complex(wp)                             :: I(4,4)

  ! Rename the parameters that describe the spin-active properties
  associate(M  => this % M_a,            &
            P  => this % polarization_a, &
            Q  => this % spinmixing_a    )

  ! Rename the Riccati parameters in the material to the left
  associate(g0   => a % g, &
            gt0  => a % gt,&
            dg0  => a % dg,&
            dgt0 => a % dgt)

  ! Calculate the spin-active terms in the interface current
  I = 0.25 * this%conductance_a * spinactive_current(g1, gt1, dg1, dgt1, g0, gt0, dg0, dgt0, M, P, Q)

  ! Calculate the deviation from the boundary condition
  r1  = r1  + (pauli0 - g1*gt1) * (I(1:2,3:4) - I(1:2,1:2)*g1)
  rt1 = rt1 + (pauli0 - gt1*g1) * (I(3:4,1:2) - I(3:4,3:4)*gt1)

  end associate
  end associate
end subroutine

pure subroutine spinactive_interface_equation_b(this, b, g2, gt2, dg2, dgt2, r2, rt2)
  ! Calculate the spin-active terms in the right boundary condition, and update the residuals.
  class(conductor), target, intent(in   ) :: this
  type(green),              intent(in   ) :: b
  type(spin),               intent(in   ) :: g2, gt2, dg2, dgt2
  type(spin),               intent(inout) :: r2, rt2

  type(spin)                              :: N2, Nt2
  type(spin)                              :: N3, Nt3
  complex(wp)                             :: L(4,4)
  complex(wp)                             :: R(4,4)
  complex(wp)                             :: I(4,4)

  ! Rename the parameters that describe the spin-active properties
  associate(M => this % M_b,            &
            P => this % polarization_b, &
            Q => this % spinmixing_b    )

  ! Rename the Riccati parameters in the material to the right
  associate(g3   => b % g, &
            gt3  => b % gt,&
            dg3  => b % dg,&
            dgt3 => b % dgt)

  ! Calculate the spin-active terms in the interface current
  I = 0.25 * this%conductance_b * spinactive_current(g2, gt2, dg2, dgt2, g3, gt3, dg3, dgt3, M, P, Q)

  ! Calculate the deviation from the boundary condition
  r2  = r2  - (pauli0 - g2*gt2) * (I(1:2,3:4) - I(1:2,1:2)*g2)
  rt2 = rt2 - (pauli0 - gt2*g2) * (I(3:4,1:2) - I(3:4,3:4)*gt2)

  end associate
  end associate
end subroutine

pure function spinactive_current(g0, gt0, dg0, dgt0, g1, gt1, dg1, dgt1, M, P, Q) result(I)
  ! Calculate the matrix current at a spin-active interface.
  type(spin),  intent(in) :: g0, gt0, dg0, dgt0
  type(spin),  intent(in) :: g1, gt1, dg1, dgt1
  real(wp),    intent(in) :: P, Q
  complex(wp), intent(in) :: M(4,4)
  complex(wp)             :: I(4,4)

  type(spin)             :: N0, Nt0
  type(spin)             :: N1, Nt1
  complex(wp)            :: GM0(4,4)
  complex(wp)            :: GM1(4,4)

  ! Calculate the normalization matrices
  N0  = spin_inv( pauli0 - g0*gt0 )
  Nt0 = spin_inv( pauli0 - gt0*g0 )
  N1  = spin_inv( pauli0 - g1*gt1 )
  Nt1 = spin_inv( pauli0 - gt1*g1 )

  ! Calculate the 4×4 Green's function in the left material
  GM0(1:2,1:2) = (+1.0_wp) * N0  * (pauli0 + g0*gt0)
  GM0(1:2,3:4) = (+2.0_wp) * N0  * g0
  GM0(3:4,1:2) = (-2.0_wp) * Nt0 * gt0
  GM0(3:4,3:4) = (-1.0_wp) * Nt0 * (pauli0 + gt0*g0)

  ! Calculate the 4×4 Green's function in the right material
  GM1(1:2,1:2) = (+1.0_wp) * N1  * (pauli0 + g1*gt1)
  GM1(1:2,3:4) = (+2.0_wp) * N1  * g1
  GM1(3:4,1:2) = (-2.0_wp) * Nt1 * gt1
  GM1(3:4,3:4) = (-1.0_wp) * Nt1 * (pauli0 + gt1*g1)

  ! Calculate the spin-active terms in the interface current
  I = spinactive_current1(GM0, GM1, M, P, Q)

contains
  pure function spinactive_current1(G0, G1, M, P, Q) result(I)
    ! Calculate the first-order matrix current. Note that we subtract the Kupriyanov-Lukichev term,
    ! so that this function can be called after the spin-inactive boundary condition without problems.
    complex(wp), intent(in) :: G0(4,4)
    complex(wp), intent(in) :: G1(4,4)
    complex(wp), intent(in) :: M(4,4)
    real(wp),    intent(in) :: P
    real(wp),    intent(in) :: Q
    complex(wp)             :: I(4,4)

    I = commutator(G0, spinactive_current1_polarization(M, G1, P) &
                     + spinactive_current1_spinmixing(M, P, Q) - G1)
  end function

  pure function spinactive_current1_polarization(M, G, P) result(F)
    ! Calculate the first-order polarization terms in the matrix current commutator.
    real(wp),    intent(in) :: P
    complex(wp), intent(in) :: G(4,4)
    complex(wp), intent(in) :: M(4,4)
    complex(wp)             :: F(4,4)

    real(wp)                :: Pr, Pp, Pm

    Pr = sqrt(1 - P**2 + eps)
    Pp = 1 + Pr
    Pm = 1 - Pr

    F  = G + (P/Pp) * anticommutator(M,G) + (Pm/Pp) * matmul(M,matmul(G,M))
  end function

  pure function spinactive_current1_spinmixing(M, P, Q) result(F)
    ! Calculate the first-order spin-mixing terms in the matrix current commutator.
    real(wp),    intent(in) :: P
    real(wp),    intent(in) :: Q
    complex(wp), intent(in) :: M(4,4)
    complex(wp)             :: F(4,4)

    real(wp)                :: Pp

    Pp = 1 + sqrt(1 - P**2 + eps)

    F = ((0,-2)*Q/Pp) * M
  end function
end function
