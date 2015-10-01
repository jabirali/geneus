! This submodule is included by conductor.f, and contains the equations which model spin-active tunneling and vacuum interfaces.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-09-01
! Updated: 2015-09-01

pure subroutine spinactive_update_prehook(this)
  ! Updates the internal variables associated with spin-active interfaces.
  class(conductor), intent(inout) :: this 
  real(dp)                        :: Pr

  if (allocated(this % magnetization_a) .and. .not. this % reflecting_a) then
    if (norm2(this % magnetization_a) < 1e-10) then
      ! Deallocate negligible magnetizations
      deallocate(this % magnetization_a)
    else
      ! Rescale the magnetization to a unit vector
      if (abs(norm2(this % magnetization_a) - 1) > 1e-10) then
        this % magnetization_a = this % magnetization_a/(norm2(this % magnetization_a) + 1e-16)
      end if

      ! Rename the relevant variables
      associate(G0 => this % conductance_a,   &
                GC => this % GC_a,            &
                GF => this % GF_a,            &
                GM => this % GM_a,            &
                H  => this % magnetization_a, & 
                P  => this % polarization_a,  &
                Q  => this % spinmixing_a,    &
                M  => this % M_a)

      ! Calculate the elements of the magnetization matrix diag(h·σ,h·σ*)
      M(1:2,1:2) = H(1) * pauli1 + H(2) * pauli2 + H(3) * pauli3
      M(3:4,3:4) = H(1) * pauli1 - H(2) * pauli2 + H(3) * pauli3

      ! Calculate the conductances associated with the spin-active interface
      if (associated(this % material_a)) then
        ! Tunneling interfaces: everything is normalized to the tunneling conductance
        Pr = sqrt(1 - P**2 + 1e-16)
        GF = 0.25 * G0 * P/(1 + Pr)
        GC = 0.25 * G0 * (1 - Pr)/(1 + Pr)
        GM = 0.25 * G0 * (-2*i*Q)/(1 + Pr)
      else
        ! Vacuum interfaces: spin-mixing term is normalized to the normal conductance
        GF = 0.00
        GC = 0.00
        GM = 0.25 * (-2*i*Q)
      end if

      end associate
    end if
  end if

  if (allocated(this % magnetization_b) .and. .not. this % reflecting_b) then
    if (norm2(this % magnetization_b) < 1e-10) then
      ! Deallocate negligible magnetizations
      deallocate(this % magnetization_b)
    else
      ! Rescale the magnetization to a unit vector
      if (abs(norm2(this % magnetization_b) - 1) > 1e-10) then
        this % magnetization_b = this % magnetization_b/(norm2(this % magnetization_b) + 1e-16)
      end if

      ! Rename the relevant variables
      associate(G0 => this % conductance_b,   &
                GC => this % GC_b,            &
                GF => this % GF_b,            &
                GM => this % GM_b,            &
                H  => this % magnetization_b, & 
                P  => this % polarization_b,  &
                Q  => this % spinmixing_b,    &
                M  => this % M_b)

      ! Calculate the elements of the magnetization matrix diag(m·σ,m·σ*)
      M(1:2,1:2) = H(1) * pauli1 + H(2) * pauli2 + H(3) * pauli3
      M(3:4,3:4) = H(1) * pauli1 - H(2) * pauli2 + H(3) * pauli3

      ! Calculate the conductances associated with the spin-active interface
      if (associated(this % material_b)) then
        ! Tunneling interfaces: everything is normalized to the tunneling conductance
        Pr = sqrt(1 - P**2 + 1e-16)
        GF = 0.25 * G0 * P/(1 + Pr)
        GC = 0.25 * G0 * (1 - Pr)/(1 + Pr)
        GM = 0.25 * G0 * (-2*i*Q)/(1 + Pr)
      else
        ! Vacuum interfaces: spin-mixing term is normalized to the normal conductance
        GF = 0.00
        GC = 0.00
        GM = 0.25 * (-2*i*Q)
      end if

      end associate
    end if
  end if
end subroutine

pure subroutine spinactive_interface_equation_a(this, a, g1, gt1, dg1, dgt1, r1, rt1)
  ! Calculate the spin-active terms in the left boundary condition, and update the residuals.
  class(conductor), target, intent(in   ) :: this
  type(green),              intent(in   ) :: a
  type(spin),               intent(in   ) :: g1, gt1, dg1, dgt1
  type(spin),               intent(inout) :: r1, rt1
  type(spin)                              :: N0, Nt0
  type(spin)                              :: N1, Nt1
  complex(dp)                             :: L(4,4)
  complex(dp)                             :: R(4,4)
  complex(dp)                             :: I(4,4)

  ! Rename the parameters that describe the spin-active properties
  associate(GC => this % GC_a,&
            GF => this % GF_a,&
            GM => this % GM_a,&
            M  => this % M_a)

  ! Rename the Riccati parameters in the material to the left
  associate(g0   => a % g, &
            gt0  => a % gt,&
            dg0  => a % dg,&
            dgt0 => a % dgt)

  ! Calculate the normalization matrices
  N0  = spin_inv( pauli0 - g0*gt0 )
  Nt0 = spin_inv( pauli0 - gt0*g0 )
  N1  = spin_inv( pauli0 - g1*gt1 )
  Nt1 = spin_inv( pauli0 - gt1*g1 )

  ! Calculate the 4×4 Green's function in the left material
  L(1:2,1:2) = (+1.0_dp) * N0  * (pauli0 + g0*gt0)
  L(1:2,3:4) = (+2.0_dp) * N0  * g0
  L(3:4,1:2) = (-2.0_dp) * Nt0 * gt0
  L(3:4,3:4) = (-1.0_dp) * Nt0 * (pauli0 + gt0*g0)

  ! Calculate the 4×4 Green's function in the right material
  R(1:2,1:2) = (+1.0_dp) * N1  * (pauli0 + g1*gt1)
  R(1:2,3:4) = (+2.0_dp) * N1  * g1
  R(3:4,1:2) = (-2.0_dp) * Nt1 * gt1
  R(3:4,3:4) = (-1.0_dp) * Nt1 * (pauli0 + gt1*g1)

  ! Calculate the spin-active terms in the interface current
  I = GC * (matmul(R,matmul(M,matmul(L,M)))    &
           -matmul(M,matmul(L,matmul(M,R))))   &
    + GF * (matmul(R,matmul(L,M)+matmul(M,L))  &
           -matmul(matmul(L,M)+matmul(M,L),R)) &
    + GM * (matmul(R,M) - matmul(M,R))

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
  complex(dp)                             :: L(4,4)
  complex(dp)                             :: R(4,4)
  complex(dp)                             :: I(4,4)

  ! Rename the parameters that describe the spin-active properties
  associate(GC => this % GC_b,&
            GF => this % GF_b,&
            GM => this % GM_b,&
            M  => this % M_b)

  ! Rename the Riccati parameters in the material to the right
  associate(g3   => b % g, &
            gt3  => b % gt,&
            dg3  => b % dg,&
            dgt3 => b % dgt)

  ! Calculate the normalization matrices
  N2  = spin_inv( pauli0 - g2*gt2 )
  Nt2 = spin_inv( pauli0 - gt2*g2 )
  N3  = spin_inv( pauli0 - g3*gt3 )
  Nt3 = spin_inv( pauli0 - gt3*g3 )

  ! Calculate the 4×4 Green's function in the left material
  L(1:2,1:2) = (+1.0_dp) * N2  * (pauli0 + g2*gt2)
  L(1:2,3:4) = (+2.0_dp) * N2  * g2
  L(3:4,1:2) = (-2.0_dp) * Nt2 * gt2
  L(3:4,3:4) = (-1.0_dp) * Nt2 * (pauli0 + gt2*g2)

  ! Calculate the 4×4 Green's function in the right material
  R(1:2,1:2) = (+1.0_dp) * N3  * (pauli0 + g3*gt3)
  R(1:2,3:4) = (+2.0_dp) * N3  * g3
  R(3:4,1:2) = (-2.0_dp) * Nt3 * gt3
  R(3:4,3:4) = (-1.0_dp) * Nt3 * (pauli0 + gt3*g3)

  ! Calculate the spin-active terms in the interface current
  I = GC * (matmul(L,matmul(M,matmul(R,M)))    &
           -matmul(M,matmul(R,matmul(M,L))))   &
    + GF * (matmul(L,matmul(R,M)+matmul(M,R))  &
           -matmul(matmul(R,M)+matmul(M,R),L)) &
    + GM * (matmul(L,M) - matmul(M,L))

  ! Calculate the deviation from the boundary condition
  r2  = r2  - (pauli0 - g2*gt2) * (I(1:2,3:4) - I(1:2,1:2)*g2)
  rt2 = rt2 - (pauli0 - gt2*g2) * (I(3:4,1:2) - I(3:4,3:4)*gt2)

  end associate
  end associate
end subroutine
