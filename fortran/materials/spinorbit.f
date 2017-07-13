!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This submodule is included by conductor.f, and contains the equations which model spin-orbit coupling in diffusive materials.

module spinorbit_m
  use :: material_m
  use :: calculus_m
  use :: matrix_m
  use :: math_m
  use :: nambu_m
  use :: spin_m
  private

  ! Public interface
  public spinorbit, spinorbit_construct

  ! Type declarations
  type :: spinorbit
    class(material), pointer    :: material  => null()                                       ! Pointer to the material modelled by this instance
    type(spin)                  :: field(3)                                                  ! Spin-orbit coupling field, i.e. SU(2) gauge field

    ! These variables are used by internal subroutines
    type(spin)                  :: Ax,  Ay,  Az,  A2                                         ! Spin-orbit coupling matrices (the components and square)
    type(spin)                  :: Axt, Ayt, Azt, A2t                                        ! Spin-orbit coupling matrices (tilde-conjugated versions)
  contains
    procedure                   :: diffusion_equation   => spinorbit_diffusion_equation      ! Defines the Usadel diffusion equation
    procedure                   :: interface_equation_a => spinorbit_interface_equation_a    ! Boundary condition at the left  interface
    procedure                   :: interface_equation_b => spinorbit_interface_equation_b    ! Boundary condition at the right interface
    procedure                   :: update_decomposition => spinorbit_update_decomposition    ! Calculate the gauge-dependent terms in the charge current decomposition
    procedure                   :: update_prehook       => spinorbit_update_prehook          ! Code to execute before calculating the propagators
    procedure                   :: update_posthook      => spinorbit_update_posthook         ! Code to execute after  calculating the propagators
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

  impure subroutine spinorbit_update_prehook(this)
    !! Updates the internal variables associated with spin-orbit coupling.
    class(spinorbit), intent(inout) :: this 

    ! Spin-orbit coupling terms in the equations for the Riccati parameter γ
    this%Ax  = this%field(1)/sqrt(this%material%thouless)
    this%Ay  = this%field(2)/sqrt(this%material%thouless)
    this%Az  = this%field(3)/sqrt(this%material%thouless)
    this%A2  = this%Ax**2 + this%Ay**2 + this%Az**2

    ! Spin-orbit coupling terms in the equations for the Riccati parameter γ~
    this%Axt = spin(conjg(this%Ax%matrix))
    this%Ayt = spin(conjg(this%Ay%matrix))
    this%Azt = spin(conjg(this%Az%matrix))
    this%A2t = spin(conjg(this%A2%matrix))
  end subroutine

  impure subroutine spinorbit_update_posthook(this)
    ! Code to execute after running the update method of a class(material) object.
    class(spinorbit), intent(inout) :: this

    ! Add gauge-covariant terms to observables
    call this % update_decomposition
  end subroutine

  pure subroutine spinorbit_diffusion_equation(this, g, gt, dg, dgt, d2g, d2gt)
    !! Calculate the spin-orbit coupling terms in the diffusion equation, and update the second derivatives of the Riccati parameters.
    class(spinorbit), intent(in)    :: this
    type(spin),       intent(in)    :: g, gt, dg, dgt
    type(spin),       intent(inout) :: d2g, d2gt
    type(spin)                      :: N, Nt

    ! Rename the spin-orbit coupling matrices
    associate(Ax => this % Ax, Axt => this % Axt, &
              Ay => this % Ay, Ayt => this % Ayt, &
              Az => this % Az, Azt => this % Azt, &
              A2 => this % A2, A2t => this % A2t  )

    ! Calculate the normalization matrices
    N   = inverse( pauli0 - g*gt )
    Nt  = inverse( pauli0 - gt*g )

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
    !! Calculate the spin-orbit coupling terms in the left boundary condition, and update the residuals.
    class(spinorbit), target, intent(in)    :: this
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
    !! Calculate the spin-orbit coupling terms in the right boundary condition, and update the residuals.
    class(spinorbit), target, intent(in)    :: this
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

  impure subroutine spinorbit_update_decomposition(this)
    !! Calculate the spin-orbit coupling terms in the charge current decomposition,
    !! i.e. determine how the SU(2) gauge field affects the relevant components.
    !!
    !! Note that the SU(2) gauge field produces a term proportional to the triple
    !! product [f(1:3)×ft(1:3)]·σ, which is an interference term that mixes the
    !! contributions from different triplet species. We therefore cannot separate
    !! the charge current exactly into contributions from specific triplets.
    !!
    !! The solution chosen here, was to divide all f(1)ft(2) and f(2)ft(1) terms
    !! equally between the contributions attributed to x- and y-triplets, etc.
    !! The advantage of this is that (i) the total triplet current is correct,
    !! and (ii) it should give a good indication of what subspecies of triplets
    !! contribute the most and least. But beware of a too literal interpretation.
    !!
    !! @TODO: The tanh(...) has to be generalized for nonequilibrium calculations.
    !!
    !! @NOTE: This function still requires debugging; it does currently not seem
    !!        to perform its function properly for the case of nanowires (Az≠0).
    !!
    class(spinorbit),         intent(inout)  :: this
    type(spin),  dimension(1:3)              :: tripleproduct
    complex(wp), dimension(0:4)              :: f, ft
    real(wp),    dimension(:,:), allocatable :: spectral
    real(wp)                                 :: prefactor
    integer                                  :: n, m

    associate(location    => this % material % location,     &
              energy      => this % material % energy,       &
              propagators => this % material % propagator,   &
              temperature => this % material % temperature,  &
              current     => this % material % decomposition,&
              Az          => this % Az                       )

      ! Allocate workspace memory
      allocate(spectral(size(energy),1:3))

      ! Iterate over the stored propagators
      do n = 1,size(location)
        do m = 1,size(energy)
          ! This factor converts from a zero-temperature to finite-temperature spectral current
          prefactor = 16 * tanh(0.8819384944310228_wp * energy(m)/temperature)

          ! Perform the singlet/triplet decomposition of the retarded propagator and its gradient
          call propagators(m,n) % decompose(f = f, ft = ft)

          ! Calculate the triple-product [f(1:3)×ft(1:3)]·σ, and partition it by triplet type
          tripleproduct(1) = (pauli2/2.0_wp) * (f(3)*ft(1) - f(1)*ft(3)) + (pauli3/2.0_wp) * (f(1)*ft(2) - f(2)*ft(1))
          tripleproduct(2) = (pauli3/2.0_wp) * (f(1)*ft(2) - f(2)*ft(1)) + (pauli1/2.0_wp) * (f(2)*ft(3) - f(3)*ft(2))
          tripleproduct(3) = (pauli1/2.0_wp) * (f(2)*ft(3) - f(3)*ft(2)) + (pauli2/2.0_wp) * (f(3)*ft(1) - f(1)*ft(3))

          ! Calculate the contribution to the spectral charge current
          spectral(m,1:3) = prefactor * re(trace(Az*tripleproduct(1:3)))
        end do

        ! Interpolate and integrate the results, and update the current vector
        current(1,n) = current(1,n) + integrate(energy, spectral(:,1), energy(1), energy(size(energy)))
        current(2,n) = current(2,n) + integrate(energy, spectral(:,2), energy(1), energy(size(energy)))
        current(3,n) = current(3,n) + integrate(energy, spectral(:,3), energy(1), energy(size(energy)))
      end do

      ! Deallocate workspace memory
      deallocate(spectral)
    end associate
  end subroutine
end module
