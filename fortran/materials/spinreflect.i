!> Author:   Jabir Ali Ouassou
!> Date:     2015-10-01
!> Category: Materials
!>
!> This submodule is included by conductor.f, and contains the equations which model spin-mixing at fully reflecting interfaces.

pure subroutine spinreflect_update_prehook(this)
  ! Updates the internal variables associated with fully reflecting spin-active interfaces.
  class(conductor), intent(inout) :: this 

  if (this % reflecting_a) then
    associate(S1 => this % S1_a, S2 => this % S2_a, Q => this % spinmixing_a)
      ! Calculate the relevant spin-mixing coefficients
      S1 = (0.00,-0.25)*sin(Q*pi)
      S2 = (0.50, 0.00)*sin(Q*pi/2)**2
    end associate
  end if

  if (this % reflecting_b) then
    associate(S1 => this % S1_b, S2 => this % S2_b, Q => this % spinmixing_b)
      ! Calculate the relevant spin-mixing coefficients
      S1 = (0.00,-0.25)*sin(Q*pi)
      S2 = (0.50, 0.00)*sin(Q*pi/2)**2
    end associate
  end if
end subroutine

pure subroutine spinreflect_interface_equation_a(this, g1, gt1, dg1, dgt1, r1, rt1)
  class(conductor), target, intent(in)    :: this
  type(spin),               intent(in)    :: g1, gt1, dg1, dgt1
  type(spin),               intent(inout) :: r1, rt1

  type(spin)                              :: N1, Nt1
  complex(wp)                             :: A(4,4)
  complex(wp)                             :: B(4,4)
  complex(wp)                             :: C(4,4)
  complex(wp)                             :: D(4,4)
  complex(wp)                             :: R(4,4)
  complex(wp)                             :: I(4,4)

  ! Rename the parameters that describe the spin-active properties
  associate(G0 => this % conductance_a, &
            S1 => this % S1_a,          &
            S2 => this % S2_a,          &
            M  => this % M_a            )

  ! Calculate the 4×4 Green's function in the left material
  R = propagator(g1, gt1)

  ! Calculate intermediate matrices used in the equation below
  A = matmul(R,matmul(M,R)) - M
  B = S1 * A
  C = S2 * matmul(A,M)
  D = S2 * matmul(M,A)

  ! Calculate the spin-active terms in the interface current
  I = matmul(matinv4(mateye4+B+C), matmul(2*matmul(R,B)+D-C, matinv4(mateye4+B+D)))/G0

  ! Calculate the deviation from the boundary condition
  r1  = r1  - (pauli0 - g1*gt1) * (I(1:2,3:4) - I(1:2,1:2)*g1)
  rt1 = rt1 - (pauli0 - gt1*g1) * (I(3:4,1:2) - I(3:4,3:4)*gt1)

  end associate
end subroutine

pure subroutine spinreflect_interface_equation_b(this, g2, gt2, dg2, dgt2, r2, rt2)
  ! Calculate the reflecting spin-mixing terms in the right boundary condition, and update the residuals.
  class(conductor), target, intent(in)    :: this
  type(spin),               intent(in)    :: g2, gt2, dg2, dgt2
  type(spin),               intent(inout) :: r2, rt2

  type(spin)                              :: N2, Nt2
  complex(wp)                             :: A(4,4)
  complex(wp)                             :: B(4,4)
  complex(wp)                             :: C(4,4)
  complex(wp)                             :: D(4,4)
  complex(wp)                             :: L(4,4)
  complex(wp)                             :: I(4,4)

  ! Rename the parameters that describe the spin-active properties
  associate(G0 => this % conductance_b, &
            S1 => this % S1_b,          &
            S2 => this % S2_b,          &
            M  => this % M_b            )

  ! Calculate the 4×4 Green's function in the left material
  L = propagator(g2, gt2)

  ! Calculate intermediate matrices used in the equation below
  A = matmul(L,matmul(M,L)) - M
  B = S1 * A
  C = S2 * matmul(A,M)
  D = S2 * matmul(M,A)

  ! Calculate the spin-active terms in the interface current
  I = matmul(matinv4(mateye4+B+C), matmul(2*matmul(L,B)+D-C, matinv4(mateye4+B+D)))/G0

  ! Calculate the deviation from the boundary condition
  r2  = r2  + (pauli0 - g2*gt2) * (I(1:2,3:4) - I(1:2,1:2)*g2)
  rt2 = rt2 + (pauli0 - gt2*g2) * (I(3:4,1:2) - I(3:4,3:4)*gt2)

  end associate
end subroutine
