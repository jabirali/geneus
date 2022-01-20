!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This submodule is included by conductor.f, and contains the equations which model spin-orbit coupling in diffusive materials.

module spinorbit_m
  use :: condmat_m
  use :: material_m
  private

  ! Public interface
  public spinorbit, spinorbit_construct

  ! Type declarations
  type :: spinorbit
    class(material), pointer   :: material  => null()                                       !! Pointer to the material modelled by this instance
    type(spin), dimension(1:3) :: field                                                     !! Spin-orbit coupling field (SU(2) gauge field)
    type(spin)                 :: Ax,  Ay,  Az,  A2                                         !! Spin-orbit coupling matrices (the components and square)
    type(spin)                 :: Axt, Ayt, Azt, A2t                                        !! Spin-orbit coupling matrices (tilde-conjugated versions)
  contains
    procedure                  :: diffusion_equation   => spinorbit_diffusion_equation      !! Diffusion equation
    procedure                  :: diffusion_equation_a => spinorbit_diffusion_equation_a    !! Boundary condition (left)
    procedure                  :: diffusion_equation_b => spinorbit_diffusion_equation_b    !! Boundary condition (right)
    procedure                  :: update_prehook       => spinorbit_update_prehook          !! Code to execute before updates
  end type

  ! Type constructors
  interface spinorbit
    module procedure spinorbit_construct
  end interface
contains
  function spinorbit_construct(parent) result(this)
    !! Constructs a spinorbit object with a given parent material.
    type(spinorbit)         :: this
    class(material), target :: parent

    ! Save a pointer to the parent object
    this % material => parent

    ! Ensure that the spin-orbit field is zero
    this % field = spin(0)
  end function

  impure recursive subroutine spinorbit_update_prehook(this)
    !! Updates the internal variables associated with spin-orbit coupling.
    class(spinorbit), intent(inout) :: this 

    ! Spin-orbit coupling terms in the equations for the Riccati parameter γ
    this % Ax  = this % material % length * this % field(1)
    this % Ay  = this % material % length * this % field(2)
    this % Az  = this % material % length * this % field(3)
    this % A2  = this % Ax**2 + this % Ay**2 + this % Az**2

    ! Spin-orbit coupling terms in the equations for the Riccati parameter γ~
    this % Axt = conjg(this % Ax)
    this % Ayt = conjg(this % Ay)
    this % Azt = conjg(this % Az)
    this % A2t = conjg(this % A2)
  end subroutine

  pure subroutine spinorbit_diffusion_equation(this, p)
    !! Calculate the spin-orbit coupling terms in the diffusion equation,
    !! and update the second derivatives of the Riccati parameters.
    class(spinorbit), intent(in)    :: this
    type(propagator), intent(inout) :: p

    associate( Ax => this % Ax,    Axt => this % Axt, &
               Ay => this % Ay,    Ayt => this % Ayt, &
               Az => this % Az,    Azt => this % Azt, &
               A2 => this % A2,    A2t => this % A2t, &
                N => p    % N,      Nt => p    % Nt,  &
                g => p    % g,      gt => p    % gt,  &
               dg => p    % dg,    dgt => p    % dgt, &
              d2g => p    % d2g,  d2gt => p    % d2gt )

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

  pure subroutine spinorbit_diffusion_equation_a(this, p, r, rt)
    !! Calculate the spin-orbit coupling terms in the left boundary condition, and update the residuals.
    class(spinorbit), target, intent(in)    :: this
    type(propagator),         intent(in)    :: p
    type(spin),               intent(inout) :: r, rt

    associate( Az  => this % Az,  Azt => this % Azt, &
                g  => p    % g,    gt => p    %  gt  )

      ! Update the residuals
      r  = r  - (0.0_wp,1.0_wp) * (Az  * g  + g  * Azt)
      rt = rt + (0.0_wp,1.0_wp) * (Azt * gt + gt * Az )

    end associate
  end subroutine

  pure subroutine spinorbit_diffusion_equation_b(this, p, r, rt)
    !! Calculate the spin-orbit coupling terms in the right boundary condition, and update the residuals.
    class(spinorbit), target, intent(in)    :: this
    type(propagator),         intent(in)    :: p
    type(spin),               intent(inout) :: r, rt

    associate( Az  => this % Az,  Azt => this % Azt, &
                g  => p    % g,    gt => p    %  gt  )

      ! Update the residuals
      r  = r  - (0.0_wp,1.0_wp) * (Az  * g  + g  * Azt)
      rt = rt + (0.0_wp,1.0_wp) * (Azt * gt + gt * Az )

    end associate
  end subroutine
end module
