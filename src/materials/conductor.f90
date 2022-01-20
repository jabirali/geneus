!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This module defines the data type 'conductor', which models the physical state of a conductor for a discretized range
!> of positions and energies.  It has two main applications: (i) it can be used as a base type for more exotic materials,
!> such as superconductors and ferromagnets; (ii) it can be used in conjunction with such materials in hybrid structures.

module conductor_m
  use :: stdio_m
  use :: condmat_m
  use :: material_m
  use :: spinorbit_m
  use :: spinactive_m
  use :: spinscattering_m
  private

  ! Type declarations
  type, public, extends(material) :: conductor
    ! Physical fields in the material
    type(spinscattering), allocatable :: spinscattering                                       !! Spin-dependent scattering
    type(spinorbit),      allocatable :: spinorbit                                            !! Spin-orbit coupling

    ! Physical fields at the interfaces
    type(spinactive),     allocatable :: spinactive_a                                         !! Spin-active interface (left)
    type(spinactive),     allocatable :: spinactive_b                                         !! Spin-active interface (right)
  contains
    ! These methods are required by the class(material) abstract interface
    procedure                 :: construct               => conductor_construct               !! Constructs the object
    procedure                 :: initialize              => conductor_initialize              !! Initializes propagators
    procedure                 :: update_prehook          => conductor_update_prehook          !! Code to execute before updates
    procedure                 :: update_posthook         => conductor_update_posthook         !! Code to execute after  updates

    ! These methods contain the equations that describe electrical conductors
    procedure                 :: diffusion_equation      => conductor_diffusion_equation      !! Diffusion equation
    procedure                 :: diffusion_equation_a    => conductor_diffusion_equation_a    !! Boundary condition (left)
    procedure                 :: diffusion_equation_b    => conductor_diffusion_equation_b    !! Boundary condition (right)

    procedure                 :: kinetic_equation        => conductor_kinetic_equation        !! Kinetic equation
    procedure                 :: kinetic_equation_a      => conductor_kinetic_equation_a      !! Boundary condition (left)
    procedure                 :: kinetic_equation_b      => conductor_kinetic_equation_b      !! Boundary condition (right)

    ! These methods define miscellaneous utility functions
    procedure                 :: conf                    => conductor_conf                    !! Configures material parameters
  end type
contains

  !--------------------------------------------------------------------------------!
  !                        IMPLEMENTATION OF CONSTRUCTORS                          !
  !--------------------------------------------------------------------------------!

  impure subroutine conductor_construct(this)
    !! Constructs a conductor object initialized to a superconducting state.
    class(conductor), intent(inout) :: this

    ! Initialize locations
    allocate(this % location(101))
    call linspace(this % location, 0 + 1e-10_wp, 1 - 1e-10_wp)

    ! Initialize energies
    allocate(this % energy(1000))
    call linspace(this % energy(   :800), 1e-6_wp, 4.00_wp)
    call linspace(this % energy(800:   ), 4.00_wp, 30.0_wp)

    ! Allocate memory for propagators
    allocate(this % propagator(size(this % energy), size(this % location)))

    ! Initialize superconductivity
    allocate(this % correlation(size(this % location)))
    this % correlation = eps

    ! Initialize magnetism
    allocate(this % magnetization(1:3, size(this % location)))
    this % magnetization = 0

    ! Allocate memory for physical observables
    allocate(this % supercurrent(0:7,size(this % location)))
    allocate(this % lossycurrent(0:7,size(this % location)))
    allocate(this % accumulation(0:7,size(this % location)))
    allocate(this % density(size(this % energy), size(this % location), 0:7))

    ! Initialize observables
    this % supercurrent = 0
    this % lossycurrent = 0
    this % accumulation = 0
    this % density      = 0

    ! Allocate boundary condition objects
    allocate(this % spinactive_a)
    allocate(this % spinactive_b)
  end subroutine

  impure subroutine conductor_initialize(this)
    !! Define the default initializer.
    class(conductor), intent(inout) :: this
    integer                         :: n, m

    ! Initialize the Riccati parameters
    do m = 1,size(this % location)
      do n = 1,size(this % energy)
        this % propagator(n,m) = propagator( cx(this % energy(n), this % scattering), this % correlation(m) )
      end do
    end do

    ! Initialize the distribution function
    do m = 1,size(this % location)
      do n = 1,size(this % energy)
        ! Finite nonequilibrium potentials
        this % propagator(n,m) % h = &
          [                                                               &
            (f(n,+1,+1) + f(n,+1,-1) + f(n,-1,+1) + f(n,-1,-1)),          &
            (f(n,+1,+1) - f(n,+1,-1) + f(n,-1,+1) - f(n,-1,-1)) * u(1),   &
            (f(n,+1,+1) - f(n,+1,-1) + f(n,-1,+1) - f(n,-1,-1)) * u(2),   &
            (f(n,+1,+1) - f(n,+1,-1) + f(n,-1,+1) - f(n,-1,-1)) * u(3),   &
            (f(n,+1,+1) + f(n,+1,-1) - f(n,-1,+1) - f(n,-1,-1)),          &
            (f(n,+1,+1) - f(n,+1,-1) - f(n,-1,+1) + f(n,-1,-1)) * u(1),   &
            (f(n,+1,+1) - f(n,+1,-1) - f(n,-1,+1) + f(n,-1,-1)) * u(2),   &
            (f(n,+1,+1) - f(n,+1,-1) - f(n,-1,+1) + f(n,-1,-1)) * u(3)    &
          ]

        ! Transverse applied potentials
        if (this % transverse) then
          this % propagator(n,m) % h(4:7) = 0
        end if

        ! Derivatives of the above
        this % propagator(n,m) % dh  = 0
        this % propagator(n,m) % d2h = 0
      end do
    end do
  contains
    pure function f(n, c, s) result(h)
      ! Fermi distribution for a given energy (n), charge parity (c=±1), and spin (s=±1).
      integer, intent(in) :: n, c, s
      real(wp)            :: h

      associate (E  => this % energy(n),           &
                 V  => this % voltage     * c,     &
                 Vs => this % spinvoltage * c * s, &
                 T  => this % temperature,         &
                 Ts => this % spintemperature * s  )
        h = tanh(0.8819384944310228_wp * (E+V+Vs)/(T+Ts))/4
      end associate
    end function

    pure function u(m) result(r)
      ! Nonequilibrium spin-projection along the m'th basis vector.
      integer, intent(in) :: m
      real(wp)            :: r

      r = this % spinaxis(m)
    end function
  end subroutine

  !--------------------------------------------------------------------------------!
  !                     IMPLEMENTATION OF CONDUCTOR EQUATIONS                      !
  !--------------------------------------------------------------------------------!

  pure subroutine conductor_diffusion_equation(this, p, e, z)
    !! Use the diffusion equation to calculate the second-derivatives 
    !! of the Riccati parameters at an energy e and position z.
    class(conductor), intent(in)    :: this
    complex(wp),      intent(in)    :: e
    real(wp),         intent(in)    :: z
    type(propagator), intent(inout) :: p

    associate(  N => p % N,     Nt => p % Nt,  &
                g => p % g,     gt => p % gt,  &
               dg => p % dg,   dgt => p % dgt, &
              d2g => p % d2g, d2gt => p % d2gt )

      ! Calculate the second-derivatives of the Riccati parameters
      d2g  = (-2.0_wp,0.0_wp)*dg*Nt*gt*dg - (0.0_wp,2.0_wp)*e*g
      d2gt = (-2.0_wp,0.0_wp)*dgt*N*g*dgt - (0.0_wp,2.0_wp)*e*gt

      ! Calculate the contribution from a spin-orbit coupling
      if (allocated(this % spinorbit)) then
        call this % spinorbit % diffusion_equation(p)
      end if

      ! Calculate the contribution from spin-dependent scattering
      if (allocated(this % spinscattering)) then
        call this % spinscattering % diffusion_equation(p)
      end if
    end associate
  end subroutine

  pure subroutine conductor_diffusion_equation_a(this, p, a, r, rt)
    !! Calculate residuals from the boundary conditions at the left interface.
    class(conductor),          intent(in)    :: this
    type(propagator),          intent(in)    :: p, a
    type(spin),                intent(inout) :: r, rt
    complex(wp), dimension(1:4,1:4)          :: I

    ! Calculate a matrix current from the propagators
    I = (+0.50_wp) * this % spinactive_a % diffusion_current(p % retarded(), a % retarded())

    ! Calculate the deviation from the boundary condition
    associate(g => p % g, gt => p % gt, dg => p % dg, dgt => p % dgt)
      r  = dg  - (pauli0 - g*gt) * (I(1:2,3:4) - I(1:2,1:2)*g)
      rt = dgt - (pauli0 - gt*g) * (I(3:4,1:2) - I(3:4,3:4)*gt)
    end associate

    ! Gauge-dependent terms in the case of spin-orbit coupling
    if (allocated(this % spinorbit)) then
      call this % spinorbit % diffusion_equation_a(p, r, rt)
    end if
  end subroutine

  pure subroutine conductor_diffusion_equation_b(this, p, b, r, rt)
    !! Calculate residuals from the boundary conditions at the right interface.
    class(conductor),          intent(in)    :: this
    type(propagator),          intent(in)    :: p, b
    type(spin),                intent(inout) :: r, rt
    complex(wp), dimension(1:4,1:4)          :: I

    ! Calculate a matrix current from the propagators
    I = (-0.50_wp) * this % spinactive_b % diffusion_current(p % retarded(), b % retarded())

    ! Calculate the deviation from the boundary condition
    associate(g => p % g, gt => p % gt, dg => p % dg, dgt => p % dgt)
      r  = dg  - (pauli0 - g*gt) * (I(1:2,3:4) - I(1:2,1:2)*g)
      rt = dgt - (pauli0 - gt*g) * (I(3:4,1:2) - I(3:4,3:4)*gt)
    end associate

    ! Gauge-dependent terms in the case of spin-orbit coupling
    if (allocated(this % spinorbit)) then
      call this % spinorbit % diffusion_equation_b(p, r, rt)
    end if
  end subroutine

  pure subroutine conductor_kinetic_equation(this, Gp, R, z)
    !! Calculate the self-energies in the kinetic equation.
    class(conductor),                intent(in)    :: this
    type(propagator),                intent(in)    :: Gp
    complex(wp), dimension(0:7,0:7), intent(inout) :: R
    real(wp),                        intent(in)    :: z
    
    ! There are no normal-metal terms
    R = 0

    ! Spin-dependent scattering terms
    if (allocated(this % spinscattering)) then
      call this % spinscattering % kinetic_equation(Gp, R)
    end if
  end subroutine

  pure subroutine conductor_kinetic_equation_a(this, Gp, Ga, Cp, Ca)
    !! Calculate proportionality matrices for the boundary conditions at the left interface.
    class(conductor),                intent(in)  :: this
    type(propagator),                intent(in)  :: Gp, Ga
    complex(wp), dimension(0:7,0:7), intent(out) :: Cp, Ca

    ! Calculate the boundary coefficients
    call this % spinactive_a % kinetic_current(Gp, Ga, Cp, Ca)

    ! Direction of the interface normal
    Cp = +Cp
    Ca = +Ca
  end subroutine

  pure subroutine conductor_kinetic_equation_b(this, Gp, Gb, Cp, Cb)
    !! Calculate proportionality matrices for the boundary conditions at the right interface.
    class(conductor),                intent(in)  :: this
    type(propagator),                intent(in)  :: Gp, Gb
    complex(wp), dimension(0:7,0:7), intent(out) :: Cp, Cb

    ! Calculate the boundary coefficients
    call this % spinactive_b % kinetic_current(Gp, Gb, Cp, Cb)

    ! Direction of the interface normal
    Cp = -Cp
    Cb = -Cb
  end subroutine

  impure recursive subroutine conductor_update_prehook(this)
    !! Code to execute before running the update method of a class(conductor) object.
    class(conductor), intent(inout) :: this
 
    ! Discard the tunneling conductance at vacuum interfaces
    if (.not. associated(this % material_a)) then
      this % spinactive_a % conductance = 1
    end if
    if (.not. associated(this % material_b)) then
      this % spinactive_b % conductance = 1
    end if

    ! Prepare variables associated with spin-orbit coupling
    if (allocated(this % spinorbit)) then
      call this % spinorbit % update_prehook
    end if

    ! Prepare variables associated with spin-active interfaces
    call this % spinactive_a % update_prehook
    call this % spinactive_b % update_prehook

    ! Modify the type string
    this % type_string = color_yellow // 'CONDUCTOR' // color_none
  end subroutine

  impure recursive subroutine conductor_update_posthook(this)
    !! Code to execute after running the update method of a class(conductor) object.
    !! In particular, this function calculates supercurrents, dissipative currents,
    !! accumulations, and density of states, and stores the results in the object.

    class(conductor), intent(inout)          :: this
    type(nambu), allocatable                 :: gauge
    real(wp),    allocatable, dimension(:,:) :: I, J, Q
    complex(wp), allocatable, dimension(:)   :: S
    integer                                  :: n, m, k

    ! Allocate memory for the workspace
    allocate(S(size(this % energy)))
    allocate(I(size(this % energy), 0:7))
    allocate(J(size(this % energy), 0:7))
    allocate(Q(size(this % energy), 0:7))

    ! Calculate the gauge contribution
    if (allocated(this % spinorbit)) then
      allocate(gauge)
      gauge = diag(+this % spinorbit % Az  % matrix,&
                   -this % spinorbit % Azt % matrix )
    end if

    ! Simplify the namespace
    associate(E => this % energy,     &
              z => this % location,   &
              G => this % propagator, &
              D => this % density     )

      ! Iterate over positions
      do n = 1,size(z)
        ! Calculate the spectral properties at this position
        do m = 1,size(E)
          S(m)     = G(m,n) % correlation()
          Q(m,:)   = G(m,n) % accumulation()
          I(m,:)   = G(m,n) % supercurrent(gauge)
          J(m,:)   = G(m,n) % lossycurrent(gauge)
          D(m,n,:) = G(m,n) % density()
        end do

        ! Superconducting correlations depend on the cutoff
        S = S/log(2*E(size(E)))

        ! Heat and spin-heat observables depend on the energy
        do k = 4,7
          Q(:,k) = E * Q(:,k)
          I(:,k) = E * I(:,k)
          J(:,k) = E * J(:,k)
        end do

        ! Integrate the spectral observables to find the total observables
        this % correlation(n) = integrate(E, S, E(1), E(size(E)))
        do k = 0,7
          this % accumulation(k,n) = integrate(E, Q(:,k), E(1), E(size(E)))
          this % supercurrent(k,n) = integrate(E, I(:,k), E(1), E(size(E)))
          this % lossycurrent(k,n) = integrate(E, J(:,k), E(1), E(size(E)))
        end do
      end do
    end associate

    ! Reset the center-of-mass phase
    if (this % phaselock) then
      associate(g => this % correlation, i => (0,1))
        g = g / exp(i*arg(mean(g)))
      end associate
    end if

    ! Deallocate workspace memory
    deallocate(S, Q, I, J)
  end subroutine



  !--------------------------------------------------------------------------------!
  !                      IMPLEMENTATION OF UTILITY METHODS                         !
  !--------------------------------------------------------------------------------!

  impure subroutine conductor_conf(this, key, val)
    !! Configure a material property based on a key-value pair.
    use :: evaluate_m

    class(conductor), intent(inout) :: this
    character(*),     intent(in)    :: key
    character(*),     intent(in)    :: val
    real(wp)                        :: tmp

    select case(key)
      case ('conductance_a')
        call evaluate(val, this % spinactive_a % conductance)
      case ('conductance_b')
        call evaluate(val, this % spinactive_b % conductance)
      case ('resistance_a')
        call evaluate(val, tmp)
        this % spinactive_a % conductance = 1/tmp
      case ('resistance_b')
        call evaluate(val, tmp)
        this % spinactive_b % conductance = 1/tmp
      case ('spinmixing_a')
        call evaluate(val, this % spinactive_a % spinmixing)
      case ('spinmixing_b')
        call evaluate(val, this % spinactive_b % spinmixing)
      case ('secondorder_a')
        call evaluate(val, this % spinactive_a % secondorder)
      case ('secondorder_b')
        call evaluate(val, this % spinactive_b % secondorder)
      case ('polarization_a')
        call evaluate(val, this % spinactive_a % polarization)
      case ('polarization_b')
        call evaluate(val, this % spinactive_b % polarization)
      case ('magnetization_a')
        call evaluate(val, this % spinactive_a % magnetization)
        this % spinactive_a % magnetization = unitvector(this % spinactive_a % magnetization)
      case ('magnetization_b')
        call evaluate(val, this % spinactive_b % magnetization)
        this % spinactive_b % magnetization = unitvector(this % spinactive_b % magnetization)
      case ('misalignment0_a')
        call evaluate(val, this % spinactive_a % misalignment0)
        this % spinactive_a % misalignment0 = unitvector(this % spinactive_a % misalignment0)
      case ('misalignment0_b')
        call evaluate(val, this % spinactive_b % misalignment0)
        this % spinactive_b % misalignment0 = unitvector(this % spinactive_b % misalignment0)
      case ('misalignment1_a')
        call evaluate(val, this % spinactive_a % misalignment1)
        this % spinactive_a % misalignment1 = unitvector(this % spinactive_a % misalignment1)
      case ('misalignment1_b')
        call evaluate(val, this % spinactive_b % misalignment1)
        this % spinactive_b % misalignment1 = unitvector(this % spinactive_b % misalignment1)
      case ('nanowire')
        call evaluate(val, tmp)
        if (.not. allocated(this % spinorbit)) then
          allocate(this % spinorbit)
          this % spinorbit = spinorbit(this)
        end if
        this % spinorbit % field(3) = this % spinorbit % field(3) + (-tmp)*pauli1
      case ('rashba')
        call evaluate(val, tmp)
        if (.not. allocated(this % spinorbit)) then
          allocate(this % spinorbit)
          this % spinorbit = spinorbit(this)
        end if
        this % spinorbit % field(1) = this % spinorbit % field(1) + (-tmp)*pauli2
        this % spinorbit % field(2) = this % spinorbit % field(2) + (+tmp)*pauli1
      case ('dresselhaus')
        call evaluate(val, tmp)
        if (.not. allocated(this % spinorbit)) then
          allocate(this % spinorbit)
          this % spinorbit = spinorbit(this)
        end if
        this % spinorbit % field(1) = this % spinorbit % field(1) + (+tmp)*pauli1
        this % spinorbit % field(2) = this % spinorbit % field(2) + (-tmp)*pauli2
      case ('gap')
        block
          real(wp), allocatable, dimension(:) :: r
          associate(g => this % correlation, z => this % location, i => (0,1))
            call evaluate(val, z, r)
            g = r * exp(i*arg(g))
          end associate
        end block
      case ('phase')
        block
          real(wp), allocatable, dimension(:) :: r
          associate(g => this % correlation, z => this % location, i => (0,1))
            call evaluate(val, z, r)
            g = abs(g) * exp(i*pi*r)
          end associate
        end block
      case ('scattering_spinflip')
        if (.not. allocated(this % spinscattering)) then
          allocate(this%spinscattering)
          this % spinscattering = spinscattering(this)
        end if
        call evaluate(val, this % spinscattering % spinflip)
      case ('scattering_spinorbit')
        if (.not. allocated(this % spinscattering)) then
          allocate(this%spinscattering)
          this % spinscattering = spinscattering(this)
        end if
        call evaluate(val, this % spinscattering % spinorbit)
      case ('depairing')
        if (.not. allocated(this % spinscattering)) then
          allocate(this%spinscattering)
          this % spinscattering = spinscattering(this)
        end if
        call evaluate(val, this % spinscattering % depairing)
      case ('zeroenergy')
        block
          logical :: tmp
          call evaluate(val, tmp)
          if (tmp) then
            deallocate(this % energy)
            allocate(this % energy(1))
            this % energy(1) = 0
          end if
        end block
      case default
        call material_conf(this, key, val)
    end select
  end subroutine
end module
