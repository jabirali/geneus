! This module defines the data type 'conductor', which models the physical state of a conductor for a discretized range
! of positions and energies.  It has two main applications: (i) it can be used as a base type for more exotic materials,
! such as superconductors and ferromagnets; (ii) it can be used in conjunction with such materials in hybrid structures.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-11
! Updated: 2015-10-04

module mod_conductor
  use mod_stdio
  use mod_math
  use mod_spin
  use mod_green
  use mod_material
  implicit none
  private

  ! Public interface
  public conductor, conductor_construct

  ! Type declarations
  type, extends(material) :: conductor
    ! These parameters control the physical characteristics of the material
    logical                   :: transparent_a           =  .false.                           ! Whether the left  interface is completely transparent
    logical                   :: transparent_b           =  .false.                           ! Whether the right interface is completely transparent
    logical                   :: reflecting_a            =  .false.                           ! Whether the left  interface is completely reflecting
    logical                   :: reflecting_b            =  .false.                           ! Whether the right interface is completely reflecting
    real(wp)                  :: conductance_a           =  0.30_wp                           ! Normalized conductance at the left  interface
    real(wp)                  :: conductance_b           =  0.30_wp                           ! Normalized conductance at the right interface
    real(wp)                  :: spinmixing_a            =  0.00_wp                           ! Normalized spin-mixing at the left  interface
    real(wp)                  :: spinmixing_b            =  0.00_wp                           ! Normalized spin-mixing at the right interface
    real(wp)                  :: polarization_a          =  0.00_wp                           ! Spin-polarization at the left  interface
    real(wp)                  :: polarization_b          =  0.00_wp                           ! Spin-polarization at the right interface

    ! These parameters represent the physical fields in the material
    real(wp)                  :: depairing               =  0.00_wp                           ! Magnetic orbital depairing
    type(spin),   allocatable :: spinorbit(:)                                                 ! Spin-orbit coupling field (spin vector)
    real(wp),     allocatable :: magnetization_a(:)                                           ! Magnetization of the left  interface (unit vector)
    real(wp),     allocatable :: magnetization_b(:)                                           ! Magnetization of the right interface (unit vector)

    ! These variables are used by internal subroutines to handle spin-active interfaces
    complex(wp),      private :: M_a(4,4)                =  0.00_wp                           ! Interface magnetization matrix in Spin-Nambu space
    complex(wp),      private :: M_b(4,4)                =  0.00_wp                           ! Interface magnetization matrix in Spin-Nambu space
    complex(wp),      private :: GM_a, GM_b                                                   ! Normalized spin-mixing G_ϕ /G_T0 at the interfaces
    real(wp),         private :: GF_a, GF_b                                                   ! Normalized spin-filter G_MR/G_T0 at the interfaces
    real(wp),         private :: GC_a, GC_b                                                   ! Normalized correction  G_T1/G_T0 at the interfaces
    complex(wp),      private :: S1_a, S1_b                                                   ! Spin-mixing prefactor proportional to sin(ϕ)
    complex(wp),      private :: S2_a, S2_b                                                   ! Spin-mixing prefactor proportional to sin(ϕ/2)²
  
    ! These variables are used by internal subroutines to handle spin-orbit coupling
    type(spin),       private :: Ax,  Ay,  Az,  A2                                            ! Spin-orbit coupling matrices (the components and square)
    type(spin),       private :: Axt, Ayt, Azt, A2t                                           ! Spin-orbit coupling matrices (tilde-conjugated versions)
  contains
    ! These methods are required by the class(material) abstract interface
    procedure                 :: init                    => conductor_init                    ! Initializes the propagators
    procedure                 :: interface_equation_a    => conductor_interface_equation_a    ! Boundary condition at the left  interface
    procedure                 :: interface_equation_b    => conductor_interface_equation_b    ! Boundary condition at the right interface
    procedure                 :: update_prehook          => conductor_update_prehook          ! Code to execute before calculating the propagators
    procedure                 :: update_posthook         => conductor_update_posthook         ! Code to execute after  calculating the propagators

    ! These methods contain the equations that describe electrical conductors
    procedure                 :: diffusion_equation      => conductor_diffusion_equation      ! Defines the Usadel diffusion equation (conductor terms)
    procedure                 :: interface_transparent_a => conductor_interface_transparent_a ! Defines the left  boundary condition  (transparent interface)
    procedure                 :: interface_transparent_b => conductor_interface_transparent_b ! Defines the right boundary condition  (transparent interface)
    procedure                 :: interface_vacuum_a      => conductor_interface_vacuum_a      ! Defines the left  boundary condition  (vacuum interface)
    procedure                 :: interface_vacuum_b      => conductor_interface_vacuum_b      ! Defines the right boundary condition  (vacuum interface)
    procedure                 :: interface_tunnel_a      => conductor_interface_tunnel_a      ! Defines the left  boundary condition  (tunnel interface)
    procedure                 :: interface_tunnel_b      => conductor_interface_tunnel_b      ! Defines the right boundary condition  (tunnel interface)

    ! These methods contain the equations that describe spin-orbit coupling
    procedure                 :: diffusion_spinorbit     => spinorbit_diffusion_equation      ! Defines the Usadel diffusion equation (spin-orbit terms)
    procedure                 :: interface_spinorbit_a   => spinorbit_interface_equation_a    ! Defines the left  boundary condition  (spin-orbit terms)
    procedure                 :: interface_spinorbit_b   => spinorbit_interface_equation_b    ! Defines the right boundary condition  (spin-orbit terms)

    ! These methods contain the equations that describe spin-active tunneling  interfaces
    procedure                 :: interface_spinactive_a  => spinactive_interface_equation_a   ! Defines the left  boundary condition (spin-active terms)
    procedure                 :: interface_spinactive_b  => spinactive_interface_equation_b   ! Defines the right boundary condition (spin-active terms)

    ! These methods contain the equations that describe spin-active reflecting interfaces
    procedure                 :: interface_spinreflect_a => spinreflect_interface_equation_a  ! Defines the left  boundary condition (spin-active terms)
    procedure                 :: interface_spinreflect_b => spinreflect_interface_equation_b  ! Defines the right boundary condition (spin-active terms)
  end type

  ! Type constructors
  interface conductor
    module procedure conductor_construct
  end interface
contains
  !--------------------------------------------------------------------------------!
  !                      INCLUDE EXTERNAL MODULE PROCEDURES                        !
  !--------------------------------------------------------------------------------!

  include 'spinorbit.i'
  include 'spinactive.i'
  include 'spinreflect.i'

  !--------------------------------------------------------------------------------!
  !                        IMPLEMENTATION OF CONSTRUCTORS                          !
  !--------------------------------------------------------------------------------!

  pure function conductor_construct(energy, gap, thouless, scattering, points) result(this)
    ! Constructs a conductor object initialized to a superconducting state.
    type(conductor)                   :: this         ! Conductor object that will be constructed
    real(wp),    intent(in)           :: energy(:)    ! Discretized energy domain that will be used
    real(wp),    intent(in), optional :: thouless     ! Thouless energy       (default: see type declaration)
    real(wp),    intent(in), optional :: scattering   ! Imaginary energy term (default: see type declaration)
    complex(wp), intent(in), optional :: gap          ! Superconducting gap   (default: see definition below)
    integer,     intent(in), optional :: points       ! Number of positions   (default: see definition below)
    integer                           :: n, m         ! Loop variables

    ! Optional argument: Thouless energy
    if (present(thouless)) then
      this%thouless = thouless
    end if

    ! Optional argument: imaginary energy
    if (present(scattering)) then
      this%scattering = scattering
    end if

    ! Allocate memory (if necessary)
    if (.not. allocated(this%propagator)) then
      if (present(points)) then
        allocate(this%propagator(size(energy), points))
        allocate(this%energy(size(energy)))
        allocate(this%location(points))
      else
        allocate(this%propagator(size(energy), 150))
        allocate(this%energy(size(energy)))
        allocate(this%location(150))
      end if
    end if

    ! Initialize energy and position arrays
    this%energy   = energy
    this%location = [ ((real(n,kind=wp)/(size(this%location)-1)), n=0,size(this%location)-1) ]

    ! Initialize the state
    if (present(gap)) then
      call this%init( gap )
    else
      call this%init( cx(1.0_wp,0.0_wp) )
    end if
  end function

  pure subroutine conductor_init(this, gap)
    ! Define the default initializer.
    class(conductor), intent(inout) :: this
    complex(wp),      intent(in   ) :: gap
    integer                         :: n, m

    do m = 1,size(this%location)
      do n = 1,size(this%energy)
        this%propagator(n,m) = green( cx(this%energy(n),this%scattering), gap )
      end do
    end do
  end subroutine

  !--------------------------------------------------------------------------------!
  !                     IMPLEMENTATION OF CONDUCTOR EQUATIONS                      !
  !--------------------------------------------------------------------------------!

  pure subroutine conductor_diffusion_equation(this, e, z, g, gt, dg, dgt, d2g, d2gt)
    ! Use the diffusion equation to calculate the second-derivatives of the Riccati parameters at energy e and point z.
    class(conductor), intent(in   ) :: this
    complex(wp),      intent(in   ) :: e
    real(wp),         intent(in   ) :: z
    type(spin),       intent(in   ) :: g, gt, dg, dgt
    type(spin),       intent(inout) :: d2g, d2gt
    type(spin)                      :: N, Nt

    ! Calculate the normalization matrices
    N   = spin_inv( pauli0 - g*gt )
    Nt  = spin_inv( pauli0 - gt*g )

    ! Calculate the second-derivatives of the Riccati parameters
    d2g  = (-2.0_wp,0.0_wp)*dg*Nt*gt*dg - (0.0_wp,2.0_wp)*e*g
    d2gt = (-2.0_wp,0.0_wp)*dgt*N*g*dgt - (0.0_wp,2.0_wp)*e*gt

    ! Calculate the contribution from a spin-orbit coupling
    if (allocated(this%spinorbit)) then
      call this%diffusion_spinorbit(g, gt, dg, dgt, d2g, d2gt)
    end if

    ! Calculate the contribution from orbital magnetic depairing
    if (this%depairing /= 0.0_wp) then
      d2g  = d2g  + (this%depairing**2/this%thouless)*(2.0_wp*N  - pauli0)*g
      d2gt = d2gt + (this%depairing**2/this%thouless)*(2.0_wp*Nt - pauli0)*gt
    end if
  end subroutine

  pure subroutine conductor_interface_equation_a(this, a, g, gt, dg, dgt, r, rt)
      ! Calculate residuals from the boundary conditions at the left interface.
      class(conductor),          intent(in   ) :: this
      type(green),               intent(in   ) :: a
      type(spin),                intent(in   ) :: g, gt, dg, dgt
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
        call this%interface_spinorbit_a(g, gt, dg, dgt, r, rt)
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
      class(conductor),          intent(in   ) :: this
      type(green),               intent(in   ) :: b
      type(spin),                intent(in   ) :: g, gt, dg, dgt
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
        call this%interface_spinorbit_b(g, gt, dg, dgt, r, rt)
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
    class(conductor), intent(in   ) :: this
    type(spin),       intent(in   ) :: g1, gt1, dg1, dgt1
    type(spin),       intent(inout) :: r1, rt1

    r1  = dg1
    rt1 = dgt1
  end subroutine

  pure subroutine conductor_interface_vacuum_b(this, g2, gt2, dg2, dgt2, r2, rt2)
    ! Defines a vacuum boundary condition for the right interface.
    class(conductor), intent(in   ) :: this
    type(spin),       intent(in   ) :: g2, gt2, dg2, dgt2
    type(spin),       intent(inout) :: r2, rt2

    r2  = dg2
    rt2 = dgt2
  end subroutine

  pure subroutine conductor_interface_transparent_a(this, a, g1, gt1, dg1, dgt1, r1, rt1)
    ! Defines a transparent boundary condition for the left interface.
    class(conductor), intent(in   ) :: this
    type(green),      intent(in   ) :: a
    type(spin),       intent(in   ) :: g1, gt1, dg1, dgt1
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
    class(conductor),          intent(in   ) :: this
    type(green),               intent(in   ) :: b
    type(spin),                intent(in   ) :: g2, gt2, dg2, dgt2
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
    class(conductor),          intent(in   ) :: this
    type(green),               intent(in   ) :: a
    type(spin),                intent(inout) :: r1, rt1
    type(spin),                intent(in   ) :: g1, gt1, dg1, dgt1
    type(spin)                               :: N0, Nt0

    ! Rename the Riccati parameters in the material to the left
    associate(g0   => a%g, &
              gt0  => a%gt,&
              dg0  => a%dg,&
              dgt0 => a%dgt)

    ! Calculate the normalization matrices
    N0  = spin_inv( pauli0 - g0*gt0 )
    Nt0 = spin_inv( pauli0 - gt0*g0 )

    ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
    r1  = dg1  - this%conductance_a*( pauli0 - g1*gt0 )*N0*(  g1  - g0  )
    rt1 = dgt1 - this%conductance_a*( pauli0 - gt1*g0 )*Nt0*( gt1 - gt0 )

    end associate
  end subroutine

  pure subroutine conductor_interface_tunnel_b(this, b, g2, gt2, dg2, dgt2, r2, rt2)
    ! Defines a tunneling boundary condition for the right interface.
    class(conductor),          intent(in   ) :: this
    type(green),               intent(in   ) :: b
    type(spin),                intent(inout) :: r2, rt2
    type(spin),                intent(in   ) :: g2, gt2, dg2, dgt2
    type(spin)                               :: N3, Nt3

    ! Rename the Riccati parameters in the material to the right
    associate(g3   => b%g, &
              gt3  => b%gt,&
              dg3  => b%dg,&
              dgt3 => b%dgt)
  
    ! Calculate the normalization matrices
    N3  = spin_inv( pauli0 - g3*gt3 )
    Nt3 = spin_inv( pauli0 - gt3*g3 )

    ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
    r2  = dg2  - this%conductance_b*( pauli0 - g2*gt3 )*N3*(  g3  - g2  )
    rt2 = dgt2 - this%conductance_b*( pauli0 - gt2*g3 )*Nt3*( gt3 - gt2 )

    end associate
  end subroutine

  impure subroutine conductor_update_prehook(this)
    ! Code to execute before running the update method of a class(conductor) object.
    class(conductor), intent(inout) :: this

    ! Prepare variables associated with spin-orbit coupling
    call spinorbit_update_prehook(this)

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
    ! Code to execute after running the update method of a class(conductor) object.
    class(conductor), intent(inout) :: this

    continue
  end subroutine
end module
