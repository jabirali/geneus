!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This module defines the data type 'conductor', which models the physical state of a conductor for a discretized range
!> of positions and energies.  It has two main applications: (i) it can be used as a base type for more exotic materials,
!> such as superconductors and ferromagnets; (ii) it can be used in conjunction with such materials in hybrid structures.
!>
!> @TODO
!>   Move the the spinactive field into separate subobjects spinactive_a and spinactive_b of a type(spinactive). This
!>   type can be moved into a module spinactive_m, together with the associated methods diffusion_spinorbit etc. We can
!>   then have separate modules spinorbit_m and spinactive_m that the current module depends on, leading to greater 
!>   encapsulation and separation. The spinreflect.i contents should also be merged into the same type and module.
!>
!> @TODO
!>   Rename spinorbit_film/spinorbit_wire/spinorbit_bulk to rashba and dresselhaus, since these are the more common
!>   names used in the litterature, and this will make it easier for other people to learn to use the code in the future.
!>   The best way to define the Rashba coupling is likely to use a vector, i.e. "rashba = [ax,ay,az]", where we can use 
!>   A=α(σ×e) to get the A-field. We can then add another option "nanowire = T" or "dimensions = 1", so spinorbit_prehook
!>   knows it should prune away the e_x and e_y components of the A-field if true. For backwards compatibility, we should
!>   perhaps interpret "rashba = az" as "rashba = [ax,ay,az]" if no brackets are detected in the value string. Think
!>   of something smart for the dresselhaus coupling as well; check e.g. Tom's thesis regarding nanowire dresselhaus.

module conductor_m
  use :: stdio_m
  use :: math_m
  use :: spin_m
  use :: propagator_m
  use :: material_m
  use :: spinorbit_m
  use :: spinscattering_m
  private

  ! Type declarations
  type, public, extends(material) :: conductor
    ! These parameters control the physical characteristics of the material
    real(wp)                         :: spinmixing_a            =  0.00_wp                    ! Normalized spin-mixing at the left  interface
    real(wp)                         :: spinmixing_b            =  0.00_wp                    ! Normalized spin-mixing at the right interface
    real(wp)                         :: polarization_a          =  0.00_wp                    ! Spin-polarization at the left  interface
    real(wp)                         :: polarization_b          =  0.00_wp                    ! Spin-polarization at the right interface
    real(wp)                         :: secondorder_a           =  0.00_wp                    ! Second-order spin-mixing at the left  interface
    real(wp)                         :: secondorder_b           =  0.00_wp                    ! Second-order spin-mixing at the right interface

    ! These parameters represent the physical fields in the material
    real(wp)                          :: depairing              =  0.00_wp                    ! Magnetic orbital depairing
    type(spinscattering), allocatable :: spinscattering                                       ! Spin-dependent scattering
    type(spinorbit),      allocatable :: spinorbit                                            ! Spin-orbit coupling
    real(wp),             allocatable :: magnetization_a(:)                                   ! Magnetization of the left  interface (unit vector) (used for transmissions)
    real(wp),             allocatable :: magnetization_b(:)                                   ! Magnetization of the right interface (unit vector) (used for transmissions)
    real(wp),             allocatable :: misalignment_a(:)                                    ! Magnetization of the left  interface (unit vector) (used for reflections)
    real(wp),             allocatable :: misalignment_b(:)                                    ! Magnetization of the right interface (unit vector) (used for reflections)

    ! These variables are used by internal subroutines to handle spin-active interfaces
    complex(wp)                      :: M_a(4,4)                =  0.00_wp                    ! Interface magnetization matrix (Spin-Nambu space) (used for transmissions through the interface)
    complex(wp)                      :: M_b(4,4)                =  0.00_wp                    ! Interface magnetization matrix (Spin-Nambu space) (used for transmissions through the interface)
    complex(wp)                      :: M0_a(4,4)               =  0.00_wp                    ! Interface magnetization matrix (Spin-Nambu space) (used for reflections on this  side of the interface)
    complex(wp)                      :: M0_b(4,4)               =  0.00_wp                    ! Interface magnetization matrix (Spin-Nambu space) (used for reflections on this  side of the interface)
    complex(wp)                      :: M1_a(4,4)               =  0.00_wp                    ! Interface magnetization matrix (Spin-Nambu space) (used for reflections on other side of the interface)
    complex(wp)                      :: M1_b(4,4)               =  0.00_wp                    ! Interface magnetization matrix (Spin-Nambu space) (used for reflections on other side of the interface)
    complex(wp),             private :: S1_a, S1_b                                            ! Spin-mixing prefactor proportional to sin(ϕ)
    complex(wp),             private :: S2_a, S2_b                                            ! Spin-mixing prefactor proportional to sin(ϕ/2)²
  contains
    ! These methods are required by the class(material) abstract interface
    procedure                 :: init                    => conductor_init                    ! Initializes the propagators
    procedure                 :: interface_equation_a    => conductor_interface_equation_a    ! Boundary condition at the left  interface
    procedure                 :: interface_equation_b    => conductor_interface_equation_b    ! Boundary condition at the right interface
    procedure                 :: update_prehook          => conductor_update_prehook          ! Code to execute before calculating the propagators
    procedure                 :: update_posthook         => conductor_update_posthook         ! Code to execute after  calculating the propagators
    procedure                 :: update_decomposition    => conductor_update_decomposition    ! Calculates the singlet/triplet decomposition

    ! These methods contain the equations that describe electrical conductors
    procedure                 :: diffusion_equation      => conductor_diffusion_equation      ! Defines the Usadel diffusion equation (conductor terms)
    procedure                 :: interface_transparent_a => conductor_interface_transparent_a ! Defines the left  boundary condition  (transparent interface)
    procedure                 :: interface_transparent_b => conductor_interface_transparent_b ! Defines the right boundary condition  (transparent interface)
    procedure                 :: interface_vacuum_a      => conductor_interface_vacuum_a      ! Defines the left  boundary condition  (vacuum interface)
    procedure                 :: interface_vacuum_b      => conductor_interface_vacuum_b      ! Defines the right boundary condition  (vacuum interface)
    procedure                 :: interface_tunnel_a      => conductor_interface_tunnel_a      ! Defines the left  boundary condition  (tunnel interface)
    procedure                 :: interface_tunnel_b      => conductor_interface_tunnel_b      ! Defines the right boundary condition  (tunnel interface)

    ! These methods contain the equations that describe spin-active tunneling  interfaces
    procedure                 :: interface_spinactive_a  => spinactive_interface_equation_a   ! Defines the left  boundary condition (spin-active terms)
    procedure                 :: interface_spinactive_b  => spinactive_interface_equation_b   ! Defines the right boundary condition (spin-active terms)

    ! These methods contain the equations that describe spin-active reflecting interfaces
    procedure                 :: interface_spinreflect_a => spinreflect_interface_equation_a  ! Defines the left  boundary condition (spin-active terms)
    procedure                 :: interface_spinreflect_b => spinreflect_interface_equation_b  ! Defines the right boundary condition (spin-active terms)

    ! These methods define miscellaneous utility functions
    procedure                 :: conf                    => conductor_conf                    ! Configures material parameters
  end type

  ! Type constructors
  interface conductor
    module procedure conductor_construct
  end interface
contains
  !--------------------------------------------------------------------------------!
  !                      INCLUDE EXTERNAL MODULE PROCEDURES                        !
  !--------------------------------------------------------------------------------!

  include 'spinactive.i'
  include 'spinreflect.i'

  !--------------------------------------------------------------------------------!
  !                        IMPLEMENTATION OF CONSTRUCTORS                          !
  !--------------------------------------------------------------------------------!

  function conductor_construct() result(this)
    ! Constructs a conductor object initialized to a superconducting state.
    type(conductor) :: this

    ! Initialize locations
    allocate(this%location(151))
    call linspace(this%location, 0 + 1e-10_wp, 1 - 1e-10_wp)

    ! Initialize energies
    allocate(this%energy(600))
    call linspace(this%energy(   :400), 1e-6_wp, 1.50_wp)
    call linspace(this%energy(400:500), 1.50_wp, 4.50_wp)
    call linspace(this%energy(500:   ), 4.50_wp, 30.0_wp)

    ! Initialize propagators
    allocate(this%propagator(size(this%energy),size(this%location)))
    call this%init( (0.0_wp,0.0_wp) )

    ! Allocate memory for physical observables
    allocate(this % correlation(size(this % location)))
    allocate(this % supercurrent(0:7,size(this % location)))
    allocate(this % lossycurrent(0:7,size(this % location)))
    allocate(this % accumulation(0:7,size(this % location)))
    allocate(this % density(size(this % energy), size(this % location), 0:7))
  end function

  pure subroutine conductor_init(this, gap)
    ! Define the default initializer.
    class(conductor),      intent(inout) :: this
    complex(wp), optional, intent(in)    :: gap
    integer                              :: n, m

    ! Initialize the Riccati parameters
    if (present(gap)) then
      do m = 1,size(this%location)
        do n = 1,size(this%energy)
            this % propagator(n,m) = propagator( cx(this%energy(n),this%scattering_inelastic), gap )
        end do
      end do
    end if

    ! Initialize the distribution function
    do m = 1,size(this%location)
      do n = 1,size(this%energy)
        this % propagator(n,m) % h = &
          [(fermi(n) + fermi(n))/2, 0.0_wp, 0.0_wp, 0.0_wp, &
           (fermi(n) - fermi(n))/2, 0.0_wp, 0.0_wp, 0.0_wp  ]
      end do
    end do
  contains
    pure function fermi(n) result(h)
      integer, intent(in) :: n
      real(wp)            :: h

      associate(E => this % energy(n), T => this % temperature)
        h = tanh(0.8819384944310228_wp * E/T)
      end associate
    end function
  end subroutine

  !--------------------------------------------------------------------------------!
  !                     IMPLEMENTATION OF CONDUCTOR EQUATIONS                      !
  !--------------------------------------------------------------------------------!

  pure subroutine conductor_diffusion_equation(this, e, z, g, gt, dg, dgt, d2g, d2gt)
    ! Use the diffusion equation to calculate the second-derivatives of the Riccati parameters at energy e and point z.
    class(conductor), intent(in)    :: this
    complex(wp),      intent(in)    :: e
    real(wp),         intent(in)    :: z
    type(spin),       intent(in)    :: g, gt, dg, dgt
    type(spin),       intent(inout) :: d2g, d2gt
    type(spin)                      :: N, Nt

    ! Calculate the normalization matrices
    N   = inverse( pauli0 - g*gt )
    Nt  = inverse( pauli0 - gt*g )

    ! Calculate the second-derivatives of the Riccati parameters
    d2g  = (-2.0_wp,0.0_wp)*dg*Nt*gt*dg - (0.0_wp,2.0_wp)*e*g
    d2gt = (-2.0_wp,0.0_wp)*dgt*N*g*dgt - (0.0_wp,2.0_wp)*e*gt

    ! Calculate the contribution from a spin-orbit coupling
    if (allocated(this%spinorbit)) then
      call this%spinorbit%diffusion_equation(g, gt, dg, dgt, d2g, d2gt)
    end if

    ! Calculate the contribution from spin-dependent scattering
    if (allocated(this%spinscattering)) then
      call this%spinscattering%diffusion_equation(g, gt, dg, dgt, d2g, d2gt)
    end if

    ! Calculate the contribution from orbital magnetic depairing
    if (this%depairing > 0) then
      d2g  = d2g  + (this%depairing/this%thouless)*(2.0_wp*N  - pauli0)*g
      d2gt = d2gt + (this%depairing/this%thouless)*(2.0_wp*Nt - pauli0)*gt
    end if
  end subroutine

  pure subroutine conductor_interface_equation_a(this, a, g, gt, dg, dgt, r, rt)
      ! Calculate residuals from the boundary conditions at the left interface.
      class(conductor),          intent(in)    :: this
      type(propagator),          intent(in)    :: a
      type(spin),                intent(in)    :: g, gt, dg, dgt
      type(spin),                intent(inout) :: r, rt

      if (associated(this%material_a)) then
        if (this%transparent_a) then
          ! Interface is transparent
          call this%interface_transparent_a(a, g, gt, dg, dgt, r, rt)
        else if (this%reflecting_a) then
          ! Interface is reflecting
          call this%interface_vacuum_a(g, gt, dg, dgt, r, rt)
        else
          ! Interface is tunneling
          call this%interface_tunnel_a(a, g, gt, dg, dgt, r, rt)
        end if
      else
        ! Interface is vacuum
        call this%interface_vacuum_a(g, gt, dg, dgt, r, rt)
      end if

      if (allocated(this%spinorbit)) then
        ! Interface has spin-orbit coupling
        call this%spinorbit%interface_equation_a(g, gt, dg, dgt, r, rt)
      end if

      if (allocated(this%magnetization_a)) then
        ! Interface is spin-active
        if (this%reflecting_a) then
          ! Must be a reflecting interface
          call this%interface_spinreflect_a(g, gt, dg, dgt, r, rt)
        else
          ! May be a tunneling interface
          call this%interface_spinactive_a(a, g, gt, dg, dgt, r, rt)
        end if
      end if
  end subroutine

  pure subroutine conductor_interface_equation_b(this, b, g, gt, dg, dgt, r, rt)
    ! Calculate residuals from the boundary conditions at the right interface.
    class(conductor),          intent(in)    :: this
    type(propagator),          intent(in)    :: b
    type(spin),                intent(in)    :: g, gt, dg, dgt
    type(spin),                intent(inout) :: r, rt

    if (associated(this%material_b)) then
      if (this%transparent_b) then
        ! Interface is transparent
        call this%interface_transparent_b(b, g, gt, dg, dgt, r, rt)
      else if (this%reflecting_b) then
        ! Interface is reflecting
        call this%interface_vacuum_b(g, gt, dg, dgt, r, rt)
      else
        ! Interface is tunneling
        call this%interface_tunnel_b(b, g, gt, dg, dgt, r, rt)
      end if
    else
      ! Interface is vacuum
      call this%interface_vacuum_b(g, gt, dg, dgt, r, rt)
    end if

    if (allocated(this%spinorbit)) then
      ! Interface has spin-orbit coupling
      call this%spinorbit%interface_equation_b(g, gt, dg, dgt, r, rt)
    end if

    if (allocated(this%magnetization_b)) then
      ! Interface is spin-active
      if (this%reflecting_b) then
        ! Must be a reflecting interface
        call this%interface_spinreflect_b(g, gt, dg, dgt, r, rt)
      else
        ! May be a tunneling interface
        call this%interface_spinactive_b(b, g, gt, dg, dgt, r, rt)
      end if
    end if
  end subroutine

  pure subroutine conductor_interface_vacuum_a(this, g1, gt1, dg1, dgt1, r1, rt1)
    ! Defines a vacuum boundary condition for the left interface.
    class(conductor), intent(in)    :: this
    type(spin),       intent(in)    :: g1, gt1, dg1, dgt1
    type(spin),       intent(inout) :: r1, rt1

    r1  = dg1
    rt1 = dgt1
  end subroutine

  pure subroutine conductor_interface_vacuum_b(this, g2, gt2, dg2, dgt2, r2, rt2)
    ! Defines a vacuum boundary condition for the right interface.
    class(conductor), intent(in)    :: this
    type(spin),       intent(in)    :: g2, gt2, dg2, dgt2
    type(spin),       intent(inout) :: r2, rt2

    r2  = dg2
    rt2 = dgt2
  end subroutine

  pure subroutine conductor_interface_transparent_a(this, a, g1, gt1, dg1, dgt1, r1, rt1)
    ! Defines a transparent boundary condition for the left interface.
    class(conductor), intent(in)    :: this
    type(propagator), intent(in)    :: a
    type(spin),       intent(in)    :: g1, gt1, dg1, dgt1
    type(spin),       intent(inout) :: r1, rt1

    ! Rename the Riccati parameters in the material to the left
    associate(g0  => a%g,&
              gt0 => a%gt)

    ! Calculate the deviation between the materials
    r1  = g1  - g0
    rt1 = gt1 - gt0

    end associate
  end subroutine

  pure subroutine conductor_interface_transparent_b(this, b, g2, gt2, dg2, dgt2, r2, rt2)
    ! Defines a transparent boundary condition for the right interface.
    class(conductor),          intent(in)    :: this
    type(propagator),          intent(in)    :: b
    type(spin),                intent(in)    :: g2, gt2, dg2, dgt2
    type(spin),                intent(inout) :: r2, rt2

    ! Rename the Riccati parameters in the material to the right
    associate(g3  => b%g,&
              gt3 => b%gt)

    ! Calculate the deviation between the materials
    r2  = g2  - g3
    rt2 = gt2 - gt3

    end associate
  end subroutine

  pure subroutine conductor_interface_tunnel_a(this, a, g1, gt1, dg1, dgt1, r1, rt1)
    ! Defines a tunneling boundary condition for the left interface.
    class(conductor),          intent(in)    :: this
    type(propagator),          intent(in)    :: a
    type(spin),                intent(inout) :: r1, rt1
    type(spin),                intent(in)    :: g1, gt1, dg1, dgt1
    type(spin)                               :: N0, Nt0

    ! Rename the Riccati parameters in the material to the left
    associate(g0   => a%g, &
              gt0  => a%gt,&
              dg0  => a%dg,&
              dgt0 => a%dgt)

    ! Calculate the normalization matrices
    N0  = inverse( pauli0 - g0*gt0 )
    Nt0 = inverse( pauli0 - gt0*g0 )

    ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
    r1  = dg1  - this%conductance_a*( pauli0 - g1*gt0 )*N0*(  g1  - g0  )
    rt1 = dgt1 - this%conductance_a*( pauli0 - gt1*g0 )*Nt0*( gt1 - gt0 )

    end associate
  end subroutine

  pure subroutine conductor_interface_tunnel_b(this, b, g2, gt2, dg2, dgt2, r2, rt2)
    ! Defines a tunneling boundary condition for the right interface.
    class(conductor),          intent(in)    :: this
    type(propagator),          intent(in)    :: b
    type(spin),                intent(inout) :: r2, rt2
    type(spin),                intent(in)    :: g2, gt2, dg2, dgt2
    type(spin)                               :: N3, Nt3

    ! Rename the Riccati parameters in the material to the right
    associate(g3   => b%g, &
              gt3  => b%gt,&
              dg3  => b%dg,&
              dgt3 => b%dgt)

    ! Calculate the normalization matrices
    N3  = inverse( pauli0 - g3*gt3 )
    Nt3 = inverse( pauli0 - gt3*g3 )

    ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
    r2  = dg2  - this%conductance_b*( pauli0 - g2*gt3 )*N3*(  g3  - g2  )
    rt2 = dgt2 - this%conductance_b*( pauli0 - gt2*g3 )*Nt3*( gt3 - gt2 )

    end associate
  end subroutine

  impure subroutine conductor_update_prehook(this)
    ! Code to execute before running the update method of a class(conductor) object.
    class(conductor), intent(inout) :: this

    ! Usually, we normalize the spin-mixing conductance and other interface parameters to the tunneling conductance. But in
    ! the case of a vacuum interface, we wish to normalize them to the normal-state conductance instead. Since the tunneling
    ! conductance is normalized to the normal conductance, we can achieve this by defining the tunneling conductance to one.
    ! Setting the polarization to zero also disables all but the spin-mixing terms in the spin-active boundary condition.
    if (.not. associated(this % material_a)) then
      this % conductance_a  = 1.0
      this % polarization_a = 0.0
    end if
    if (.not. associated(this % material_b)) then
      this % conductance_b  = 1.0
      this % polarization_b = 0.0
    end if

    ! Prepare variables associated with spin-orbit coupling
    if (allocated(this%spinorbit)) then
      call this%spinorbit%update_prehook
    end if

    ! Prepare variables associated with spin-active tunneling  interfaces
    call spinactive_update_prehook(this)

    ! Prepare variables associated with spin-active reflecting interfaces
    call spinreflect_update_prehook(this)

    ! Modify the type string
    this%type_string = color_yellow // 'CONDUCTOR' // color_none
    if (allocated(this%spinorbit))       this%type_string = trim(this%type_string) // color_cyan   // ' [SOC]' // color_none
    if (allocated(this%magnetization_a)) this%type_string = trim(this%type_string) // color_purple // ' [SAL]' // color_none
    if (allocated(this%magnetization_b)) this%type_string = trim(this%type_string) // color_purple // ' [SAR]' // color_none
  end subroutine

  impure subroutine conductor_update_posthook(this)
    !! Code to execute after running the update method of a class(conductor) object.
    !! In particular, this function calculates supercurrents, dissipative currents,
    !! accumulations, and density of states, and stores the results in the object.
    use :: nambu_m

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
        S = S/acosh(E(size(E)))

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

    ! Deallocate workspace memory
    deallocate(S, Q, I, J)

    ! Calculate the current decomposition
    call this%update_decomposition

    ! Call the spinorbit posthook
    if (allocated(this%spinorbit)) then
      call this%spinorbit%update_posthook
    end if
  end subroutine

  impure subroutine conductor_update_decomposition(this)
    !! Calculate the singlet/triplet decomposition of the charge current in the material.
    !! @TODO: The tanh(...) has to be generalized for future nonequilibrium calculations.
    class(conductor), intent(inout) :: this
    complex(wp)                     :: f(0:3), df(0:3), ft(0:3), dft(0:3)
    real(wp),         allocatable   :: spectral(:,:)
    real(wp)                        :: prefactor(0:3)
    integer                         :: n, m

    if (size(this%energy) > 1) then
      ! Allocate memory if necessary
      if (.not. allocated(this%decomposition)) then
        allocate(this%decomposition(0:3,size(this%location)))
      end if

      ! Allocate workspace memory
      allocate(spectral(size(this%energy),0:3))

      ! Iterate over the stored propagators
      do n = 1,size(this%location)
        do m = 1,size(this%energy)
          ! This factor converts from a zero-temperature to finite-temperature spectral current, and
          ! adds the signs that distinguish between the (positive) singlet and (negative) triplet current
          prefactor = [+1,-1,-1,-1] * 8 * tanh(0.8819384944310228_wp * this%energy(m)/this%temperature)

          ! Perform the singlet/triplet decomposition of the retarded propagator and its gradient
          call this % propagator(m,n) % decompose(f = f, ft = ft, df = df, dft = dft)

          ! Calculate the contribution to the spectral charge current
          spectral(m,:) = prefactor * re(f*dft - ft*df)
        end do

        ! Interpolate and integrate the results, and update the current vector
        this%decomposition(0,n) = integrate(this%energy, spectral(:,0), this%energy(1), this%energy(ubound(this%energy,1)))
        this%decomposition(1,n) = integrate(this%energy, spectral(:,1), this%energy(1), this%energy(ubound(this%energy,1)))
        this%decomposition(2,n) = integrate(this%energy, spectral(:,2), this%energy(1), this%energy(ubound(this%energy,1)))
        this%decomposition(3,n) = integrate(this%energy, spectral(:,3), this%energy(1), this%energy(ubound(this%energy,1)))
      end do

      ! Deallocate workspace memory
      deallocate(spectral)
    end if
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
      case ('transparent_a')
        call evaluate(val, this % transparent_a)
      case ('transparent_b')
        call evaluate(val, this % transparent_b)
      case ('reflecting_a')
        call evaluate(val, this % reflecting_a)
      case ('reflecting_b')
        call evaluate(val, this % reflecting_b)
      case ('conductance_a')
        call evaluate(val, this % conductance_a)
      case ('conductance_b')
        call evaluate(val, this % conductance_b)
      case ('resistance_a')
        call evaluate(val, tmp)
        this % conductance_a = 1/tmp
      case ('resistance_b')
        call evaluate(val, tmp)
        this % conductance_b = 1/tmp
      case ('spinmixing_a')
        call evaluate(val, this % spinmixing_a)
      case ('spinmixing_b')
        call evaluate(val, this % spinmixing_b)
      case ('secondorder_a')
        call evaluate(val, this % secondorder_a)
      case ('secondorder_b')
        call evaluate(val, this % secondorder_b)
      case ('polarization_a')
        call evaluate(val, this % polarization_a)
      case ('polarization_b')
        call evaluate(val, this % polarization_b)
      case ('magnetization_a')
        if (.not. allocated(this % magnetization_a)) then
          allocate(this % magnetization_a(3))
        end if
        call evaluate(val, this % magnetization_a)
        if (norm2(this % magnetization_a) > sqrt(eps)) then
          this % magnetization_a = unitvector(this % magnetization_a)
        else
          deallocate(this % magnetization_a)
        end if
      case ('magnetization_b')
        if (.not. allocated(this % magnetization_b)) then
          allocate(this % magnetization_b(3))
        end if
        call evaluate(val, this % magnetization_b)
        if (norm2(this % magnetization_b) > sqrt(eps)) then
          this % magnetization_b = unitvector(this % magnetization_b)
        else
          deallocate(this % magnetization_b)
        end if
      case ('misalignment_a')
        if (.not. allocated(this % misalignment_a)) then
          allocate(this % misalignment_a(3))
        end if
        call evaluate(val, this % misalignment_a)
        if (norm2(this % misalignment_a) > sqrt(eps)) then
          this % misalignment_a = unitvector(this % misalignment_a)
        else
          deallocate(this % misalignment_a)
        end if
      case ('misalignment_b')
        if (.not. allocated(this % misalignment_b)) then
          allocate(this % misalignment_b(3))
        end if
        call evaluate(val, this % misalignment_b)
        if (norm2(this % misalignment_b) > sqrt(eps)) then
          this % misalignment_b = unitvector(this % misalignment_b)
        else
          deallocate(this % misalignment_b)
        end if
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
      case ('depairing')
        call evaluate(val, this % depairing)
      case ('gap')
        block
          integer  :: index
          real(wp) :: gap, phase
          index = scan(val,',')
          if (index > 0) then
            call evaluate(val(1:index-1), gap)
            call evaluate(val(index+1: ), phase)
          else
            call evaluate(val, gap)
            phase = 0
          end if
          call this % init( gap = gap*exp((0.0,1.0)*pi*phase) )
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
