!> Author:   Jabir Ali Ouassou
!> Date:     2015-10-01
!> Category: Materials
!>
!> This submodule is included by conductor.f, and contains the equations which model spin-orbit coupling in diffusive materials.

module spinorbit_m
  use :: material_m
  use :: calculus_m
  use :: matrix_m
  use :: math_m
  use :: spin_m
  private

  ! Public interface
  public spinorbit, spinorbit_construct

  ! Type declarations
  type :: spinorbit
    class(material), pointer    :: material  => null()                                       ! Pointer to the material modelled by this instance
    type(spin)                  :: field(3)                                                  ! Spin-orbit coupling field, i.e. SU(2) gauge field

    ! These variables are used by internal subroutines
    complex(wp), dimension(4,4) :: sigma0, sigma1, sigma2, sigma3                            ! Pauli matrices that are used internally in this object
    type(spin)                  :: Ax,  Ay,  Az,  A2                                         ! Spin-orbit coupling matrices (the components and square)
    type(spin)                  :: Axt, Ayt, Azt, A2t                                        ! Spin-orbit coupling matrices (tilde-conjugated versions)
  contains
    procedure                   :: diffusion_equation   => spinorbit_diffusion_equation      ! Defines the Usadel diffusion equation
    procedure                   :: interface_equation_a => spinorbit_interface_equation_a    ! Boundary condition at the left  interface
    procedure                   :: interface_equation_b => spinorbit_interface_equation_b    ! Boundary condition at the right interface
    procedure                   :: update_current       => spinorbit_update_current          ! Calculate the gauge-dependent terms in the current
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

    ! Construct the necessary 4×4 basis matrices
    this % sigma0 = diag(+pauli0%matrix, -pauli0%matrix)
    this % sigma1 = diag(+pauli1%matrix, -pauli1%matrix)
    this % sigma2 = diag(+pauli2%matrix, +pauli2%matrix)
    this % sigma3 = diag(+pauli3%matrix, -pauli3%matrix)
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

    call this % update_current
  end subroutine

  pure subroutine spinorbit_diffusion_equation(this, g, gt, dg, dgt, d2g, d2gt)
    !! Calculate the spin-orbit coupling terms in the diffusion equation, and update the second derivatives of the Riccati parameters.
    class(spinorbit), intent(in)    :: this
    type(spin),       intent(in)    :: g, gt, dg, dgt
    type(spin),       intent(inout) :: d2g, d2gt
    type(spin)                      :: N,  Nt

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

  impure subroutine spinorbit_update_current(this)
    !! Calculate the spin-orbit coupling terms in the charge and spin currents, 
    !! i.e. determine how the SU(2) gauge field affects the relevant currents.
    class(spinorbit), intent(inout)  :: this
    real(wp),         allocatable    :: spectral(:,:)
    real(wp)                         :: prefactor
    complex(wp),      dimension(4,4) :: G, Gt, A, K
    integer                          :: n, m

    associate(location    => this % material % location,    &
              energy      => this % material % energy,      &
              propagators => this % material % propagator,  &
              temperature => this % material % temperature, &
              current     => this % material % current      )

      ! Allocate workspace memory
      allocate(spectral(size(energy),0:3))

      ! Construct the 4×4 spin-orbit matrix
      A  = diag(+this%Az%matrix, -this%Azt%matrix)

      ! Iterate over the stored propagators
      do n = 1,size(location)
        do m = 1,size(energy)
          ! This factor converts from a zero-temperature to finite-temperature spectral current
          prefactor = 2 * tanh(0.8819384944310228_wp * energy(m)/temperature)

          ! Construct the 4×4 propagator matrices at this position and energy
          ! (The `block` prevents a segfault when compiling with IFort 16 and full 
          !  optimization; I have no idea why. It has no effect under GFortran 6.)
          block
            G  = propagators(m,n) % matrix()
            Gt = conjg(propagators(m,n) % matrixt())
          end block

          ! Calculate the corresponding spin-orbit contribution to the 4×4 spectral matrix current
          K = matmul(G,matmul(A,G)) - matmul(Gt,matmul(A,Gt))

          ! Calculate the contribution to the spectral charge and spin currents at this position
          spectral(m,0) = prefactor * im(trace(matmul(this%sigma0,K)))
          spectral(m,1) = prefactor * im(trace(matmul(this%sigma1,K)))
          spectral(m,2) = prefactor * im(trace(matmul(this%sigma2,K)))
          spectral(m,3) = prefactor * im(trace(matmul(this%sigma3,K)))
        end do

        ! Interpolate and integrate the results, and update the current vector
        current(0,n) = current(0,n) + integrate(energy, spectral(:,0), energy(1), energy(size(energy)))
        current(1,n) = current(1,n) + integrate(energy, spectral(:,1), energy(1), energy(size(energy)))
        current(2,n) = current(2,n) + integrate(energy, spectral(:,2), energy(1), energy(size(energy)))
        current(3,n) = current(3,n) + integrate(energy, spectral(:,3), energy(1), energy(size(energy)))
      end do

      ! Deallocate workspace memory
      deallocate(spectral)
    end associate
  end subroutine
end module
