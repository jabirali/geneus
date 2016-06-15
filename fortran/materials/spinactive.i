!> Author:   Jabir Ali Ouassou
!> Date:     2015-10-01
!> Category: Materials
!>
!> This submodule is included by conductor.f, and contains the equations which model spin-active tunneling and vacuum interfaces.

pure subroutine spinactive_update_prehook(this)
  !! Updates the internal variables associated with spin-active interfaces.
  class(conductor), intent(inout) :: this 

  ! Process transmission properties (both sides of interfaces)
  call update_magnetization(this % M_a, this % magnetization_a)
  call update_magnetization(this % M_b, this % magnetization_b)

  ! Default reflection properties match transmission properties
  this % M0_a = this % M_a
  this % M0_b = this % M_b
  this % M1_a = this % M_a
  this % M1_b = this % M_b

  ! Process reflection properties (this side of interfaces)
  call update_magnetization(this % M0_a, this % misalignment_a)
  call update_magnetization(this % M0_b, this % misalignment_b)

  ! Process reflection properties (other side of interfaces)
  if (associated(this % material_a)) then
    select type (other => this % material_a)
      class is (conductor)
        call update_magnetization(this % M1_a, other % misalignment_b)
    end select
  end if

  if (associated(this % material_b)) then
    select type (other => this % material_b)
      class is (conductor)
        call update_magnetization(this % M1_b, other % misalignment_a)
    end select
  end if
contains
  pure subroutine update_magnetization(matrix, vector)
    !! Updates a magnetization matrix based on the content of an allocatable magnetization vector. 
    !! If the magnetization vector is not allocated, then the magnetization matrix is not updated.
    real(wp), allocatable, intent(in)    :: vector(:)
    complex(wp),           intent(inout) :: matrix(4,4)

    if (allocated(vector)) then
      matrix(1:2,1:2) = vector(1) * pauli1 + vector(2) * pauli2 + vector(3) * pauli3
      matrix(3:4,3:4) = vector(1) * pauli1 - vector(2) * pauli2 + vector(3) * pauli3
    end if
  end subroutine
end subroutine

pure subroutine spinactive_interface_equation_a(this, a, g1, gt1, dg1, dgt1, r1, rt1)
  !! Calculate the spin-active terms in the left boundary condition, and update the residuals.
  class(conductor), target,   intent(in)    :: this
  type(propagator),           intent(in)    :: a
  type(spin),                 intent(in)    :: g1, gt1, dg1, dgt1
  type(spin),                 intent(inout) :: r1, rt1
  complex(wp), dimension(4,4)               :: GM0, GM1, I

  ! Calculate the 4×4 matrix propagators
  associate(g0 => a % g, gt0 => a % gt)
    GM0 = propagator(g0, gt0)
    GM1 = propagator(g1, gt1)
  end associate

  ! Calculate the 4×4 matrix current
  I = 0.25 * this%conductance_a &
    * spinactive_current(GM1, GM0, this%M_a, this%M0_a, this%M1_a, this%polarization_a, this%spinmixing_a, this%secondorder_a)

  ! Calculate the deviation from the boundary condition
  r1  = r1  + (pauli0 - g1*gt1) * (I(1:2,3:4) - I(1:2,1:2)*g1)
  rt1 = rt1 + (pauli0 - gt1*g1) * (I(3:4,1:2) - I(3:4,3:4)*gt1)
end subroutine

pure subroutine spinactive_interface_equation_b(this, b, g2, gt2, dg2, dgt2, r2, rt2)
  !! Calculate the spin-active terms in the right boundary condition, and update the residuals.
  class(conductor), target,   intent(in)    :: this
  type(propagator),           intent(in)    :: b
  type(spin),                 intent(in)    :: g2, gt2, dg2, dgt2
  type(spin),                 intent(inout) :: r2, rt2
  complex(wp), dimension(4,4)               :: GM2, GM3, I

  ! Calculate the 4×4 matrix propagators
  associate(g3 => b % g, gt3 => b % gt)
    GM2 = propagator(g2, gt2)
    GM3 = propagator(g3, gt3)
  end associate

  ! Calculate the 4×4 matrix current
  I = 0.25 * this%conductance_b &
    * spinactive_current(GM2, GM3, this%M_b, this%M0_b, this%M1_b, this%polarization_b, this%spinmixing_b, this%secondorder_b)

  ! Calculate the deviation from the boundary condition
  r2  = r2  - (pauli0 - g2*gt2) * (I(1:2,3:4) - I(1:2,1:2)*g2)
  rt2 = rt2 - (pauli0 - gt2*g2) * (I(3:4,1:2) - I(3:4,3:4)*gt2)
end subroutine

pure function spinactive_current(G0, G1, M, M0, M1, P, Q, R) result(I)
  !! Calculate the matrix current at an interface with spin-active properties. The equations
  !! implemented here should be valid for an arbitrary interface polarization, and up to 2nd
  !! order in the transmission probabilities and spin-mixing angles of the interface. 
  use :: matrix_m

  complex(wp), dimension(4,4), intent(in) :: G0, G1      !! Propagator matrices
  complex(wp), dimension(4,4), intent(in) :: M0, M1, M   !! Magnetization matrices 
  real(wp),                    intent(in) :: P,  Q,  R   !! Interface parameters
  complex(wp), dimension(4,4)             :: S0, S1      !! Matrix expressions
  complex(wp), dimension(4,4)             :: I           !! Matrix current

  ! Calculate the 1st-order contributions to the matrix current.  Note that we manually
  ! subtract the Kupriyanov-Lukichev term from the result, so that this function can be
  ! called after the spin-inactive boundary conditions without double-counting terms.

  ! Evaluate the 1st-order matrix functions
  S0 = spinactive_current1_transmission(G1)
  S1 = spinactive_current1_reflection()

  ! Evaluate the 1st-order matrix current
  I  = commutator(G0, S0 + S1 - G1)

  ! Calculate the 2nd-order contributions to the matrix current. Note that we make a
  ! number of simplifications in this implementation. In particular, we assume that 
  ! all interface parameters except the magnetization directions are equal on both
  ! sides of the interface. We also assume that the spin-mixing angles and tunneling
  ! probabilities of different channels have standard deviations that are much smaller
  ! than their mean values, which reduces the number of new fitting parameters to one.

  if (abs(R) > 0) then
    ! Evaluate the 1st-order matrix functions
    S1 = spinactive_current1_transmission(matmul(G1,matmul(M1,G1)) - M1) 

    ! Evaluate the 2nd-order matrix current
    I = I                                     &
      + spinactive_current2_transmission()    &
      + spinactive_current2_crossterms()      &
      + spinactive_current2_reflection()
  end if
contains
  pure function spinactive_current1_transmission(G) result(F)
    !! Calculate the 1st-order transmission terms in the matrix current commutator.
    complex(wp), dimension(4,4), intent(in) :: G
    complex(wp), dimension(4,4)             :: F

    real(wp) :: Pr, Pp, Pm
    Pr = sqrt(1 - P**2)
    Pp = 1 + Pr
    Pm = 1 - Pr

    F  = G + (P/Pp) * anticommutator(M,G) + (Pm/Pp) * matmul(M,matmul(G,M))
  end function

  pure function spinactive_current1_reflection() result(F)
    !! Calculate the 1st-order spin-mixing terms in the matrix current commutator.
    complex(wp), dimension(4,4) :: F

    F = ((0,-1)*Q) * M0
  end function

  pure function spinactive_current2_transmission() result(I)
    !! Calculate the 2nd-order transmission terms in the matrix current.
    complex(wp), dimension(4,4) :: I

    I = (0.50*R/Q) * matmul(S0,matmul(G0,S0))
  end function

  pure function spinactive_current2_reflection() result(I)
    !! Calculate the 2nd-order spin-mixing terms in the matrix current.
    complex(wp), dimension(4,4) :: I

    I = (0.25*R*Q) * commutator(G0, matmul(M0,matmul(G0,M0)))
  end function

  pure function spinactive_current2_crossterms() result(I)
    !! Calculate the 2nd-order cross-terms in the matrix current.
    complex(wp), dimension(4,4) :: I

    I = ((0.00,0.25)*R) * commutator(G0, matmul(S0,matmul(G0,M0)) + matmul(M0,matmul(G0,S0)) + S1)
  end function
end function
