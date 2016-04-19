! This submodule is included by conductor.f, and contains the equations which model spin-orbit coupling in diffusive materials.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-10-01
! Updated: 2015-10-04

pure subroutine spinorbit_update_prehook(this)
  ! Updates the internal variables associated with spin-orbit coupling.
  class(conductor), intent(inout) :: this 

  if (allocated(this%spinorbit)) then
    if (sum(this%spinorbit%norm()) < sqrt(eps)) then
      ! Negligible spin-orbit coupling
      deallocate(this%spinorbit)
    else
      ! Spin-orbit coupling terms in the equations for the Riccati parameter γ
      this%Ax  = this%spinorbit(1)/sqrt(this%thouless)
      this%Ay  = this%spinorbit(2)/sqrt(this%thouless)
      this%Az  = this%spinorbit(3)/sqrt(this%thouless)
      this%A2  = this%Ax**2 + this%Ay**2 + this%Az**2

      ! Spin-orbit coupling terms in the equations for the Riccati parameter γ~
      this%Axt = spin(conjg(this%Ax%matrix))
      this%Ayt = spin(conjg(this%Ay%matrix))
      this%Azt = spin(conjg(this%Az%matrix))
      this%A2t = spin(conjg(this%A2%matrix))
    end if
  end if
end subroutine

pure subroutine spinorbit_diffusion_equation(this, g, gt, dg, dgt, d2g, d2gt)
  ! Calculate the spin-orbit coupling terms in the diffusion equation, and update the second derivatives of the Riccati parameters.
  class(conductor), intent(in)    :: this
  type(spin),       intent(in)    :: g, gt, dg, dgt
  type(spin),       intent(inout) :: d2g, d2gt
  type(spin)                      :: N,  Nt

  ! Rename the spin-orbit coupling matrices
  associate(Ax => this % Ax, Axt => this % Axt,&
            Ay => this % Ay, Ayt => this % Ayt,&
            Az => this % Az, Azt => this % Azt,&
            A2 => this % A2, A2t => this % A2t)

  ! Calculate the normalization matrices
  N   = spin_inv( pauli0 - g*gt )
  Nt  = spin_inv( pauli0 - gt*g )

  ! Update the second derivatives of the Riccati parameters
  d2g  = d2g             + (A2 * g - g * A2t)                             &
       + (2.0_wp,0.0_wp) * (Ax * g + g * Axt) * Nt * (Axt + gt * Ax * g)  &
       + (2.0_wp,0.0_wp) * (Ay * g + g * Ayt) * Nt * (Ayt + gt * Ay * g)  &
       + (2.0_wp,0.0_wp) * (Az * g + g * Azt) * Nt * (Azt + gt * Az * g)  &
       + (0.0_wp,2.0_wp) * (Az + g * Azt * gt) * N * dg                   &
       + (0.0_wp,2.0_wp) * dg * Nt * (gt * Az * g + Azt)

  d2gt = d2gt            + (A2t * gt - gt * A2)                           &
       + (2.0_wp,0.0_wp) * (Axt * gt + gt * Ax) * N * (Ax + g * Axt * gt) &
       + (2.0_wp,0.0_wp) * (Ayt * gt + gt * Ay) * N * (Ay + g * Ayt * gt) &
       + (2.0_wp,0.0_wp) * (Azt * gt + gt * Az) * N * (Az + g * Azt * gt) &
       - (0.0_wp,2.0_wp) * (Azt + gt * Az * g) * Nt * dgt                 &
       - (0.0_wp,2.0_wp) * dgt * N * (g * Azt * gt + Az)

  end associate
end subroutine

pure subroutine spinorbit_interface_equation_a(this, g1, gt1, dg1, dgt1, r1, rt1)
  ! Calculate the spin-orbit coupling terms in the left boundary condition, and update the residuals.
  class(conductor), target, intent(in)    :: this
  type(spin),               intent(in)    :: g1, gt1, dg1, dgt1
  type(spin),               intent(inout) :: r1, rt1

  ! Rename the spin-orbit coupling matrices
  associate(Az  => this % Az,&
            Azt => this % Azt)

  ! Update the residuals
  r1  = r1  - (0.0_wp,1.0_wp) * (Az  * g1  + g1  * Azt)
  rt1 = rt1 + (0.0_wp,1.0_wp) * (Azt * gt1 + gt1 * Az )

  end associate
end subroutine

pure subroutine spinorbit_interface_equation_b(this, g2, gt2, dg2, dgt2, r2, rt2)
  ! Calculate the spin-orbit coupling terms in the right boundary condition, and update the residuals.
  class(conductor), target, intent(in)    :: this
  type(spin),               intent(in)    :: g2, gt2, dg2, dgt2
  type(spin),               intent(inout) :: r2, rt2

  ! Rename the spin-orbit coupling matrices
  associate(Az   => this % Az,&
            Azt  => this % Azt)

  ! Update the residuals
  r2  = r2  - (0.0_wp,1.0_wp) * (Az  * g2  + g2  * Azt)
  rt2 = rt2 + (0.0_wp,1.0_wp) * (Azt * gt2 + gt2 * Az )  

  end associate
end subroutine
