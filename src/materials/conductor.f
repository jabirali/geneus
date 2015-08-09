! This module defines the data type 'conductor', which models the physical state of a conductor for a discretized range
! of positions and energies.  It has two main applications: (i) it can be used as a base type for more exotic materials,
! such as superconductors and ferromagnets; (ii) it can be used in conjunction with such materials in hybrid structures. 
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-11
! Updated: 2015-08-09

module mod_conductor
  use mod_material
  implicit none

  ! Type declarations
  type, extends(material) :: conductor
    ! These parameters control the physical characteristics of the material 
    logical                   :: transparent_a           =  .false.                           ! Whether the left  interface is transparent or not
    logical                   :: transparent_b           =  .false.                           ! Whether the right interface is transparent or not
    real(dp)                  :: conductance_a           =  0.00_dp                           ! Normalized conductance of the left interface  (relative to the bulk conductance of this material)
    real(dp)                  :: conductance_b           =  0.00_dp                           ! Normalized conductance of the right interface (relative to the bulk conductance of this material)
    real(dp)                  :: polarization_a          =  0.00_dp                           ! Spin-polarization of the left interface  (range: [-1,+1])
    real(dp)                  :: polarization_b          =  0.00_dp                           ! Spin-polarization of the right interface (range: [-1,+1])
    real(dp)                  :: phaseshift_a            =  0.00_dp                           ! Spin-dependent phase shifts at the left interface
    real(dp)                  :: phaseshift_b            =  0.00_dp                           ! Spin-dependent phase shifts at the right interface

    ! These parameters represent the physical fields in the material
    real(dp),     allocatable :: magnetization_a(:)                                           ! Magnetization of the left interface  (unit vector)
    real(dp),     allocatable :: magnetization_b(:)                                           ! Magnetization of the right interface (unit vector)
    type(spin),   allocatable :: spinorbit(:)                                                 ! Spin-orbit coupling field (spin vector)

    ! These variables are used by internal subroutines to handle spin-active interfaces
    complex(dp),      private :: M_a(4,4)                =  0.00_dp                           ! Left interface magnetization matrix in Spin-Nambu space
    complex(dp),      private :: M_b(4,4)                =  0.00_dp                           ! Right interface magnetization matrix in Spin-Nambu space
    complex(dp),      private :: GP_a, GP_b                                                   ! Normalized phaseshift factors G_ϕ/G at the left and right interfaces
    real(dp),         private :: GM_a, GM_b                                                   ! Normalized magnetoconductances G_MR/G at the left and right interfaces
    real(dp),         private :: G1_a, G1_b                                                   ! Normalized correction factors G_T1/G at the left and right interfaces
  
    ! These variables are used by internal subroutines to handle spin-orbit coupling
    type(spin),       private :: Ax,  Ay,  Az,  A2                                            ! Spin-orbit coupling matrices (the components and square)
    type(spin),       private :: Axt, Ayt, Azt, A2t                                           ! Spin-orbit coupling matrices (tilde-conjugated versions)
  contains
    ! These methods are required by the class(material) abstract interface
    procedure                 :: init                    => conductor_init                    ! Initializes the Green's functions
    procedure                 :: interface_equation_a    => conductor_interface_equation_a    ! Boundary condition at the left  interface
    procedure                 :: interface_equation_b    => conductor_interface_equation_b    ! Boundary condition at the right interface
    procedure                 :: update_prehook          => conductor_update_prehook          ! Code to execute before calculating the Green's functions
    procedure                 :: update_posthook         => conductor_update_posthook         ! Code to execute after  calculating the Green's functions

    ! These methods contain the equations that describe electrical conductors
    procedure                 :: diffusion_equation      => conductor_diffusion_equation      ! Defines the Usadel diffusion equation
    procedure                 :: interface_transparent_a => conductor_interface_transparent_a ! Defines the left  boundary condition (tunnel interface)
    procedure                 :: interface_transparent_b => conductor_interface_transparent_b ! Defines the right boundary condition (tunnel interface)
    procedure                 :: interface_vacuum_a      => conductor_interface_vacuum_a      ! Defines the left  boundary condition (vacuum interface)
    procedure                 :: interface_vacuum_b      => conductor_interface_vacuum_b      ! Defines the right boundary condition (vacuum interface)
    procedure                 :: interface_tunnel_a      => conductor_interface_tunnel_a      ! Defines the left  boundary condition (tunnel interface)
    procedure                 :: interface_tunnel_b      => conductor_interface_tunnel_b      ! Defines the right boundary condition (tunnel interface)

    ! These methods contain the equations that describe spin-orbit coupling
    procedure                 :: diffusion_spinorbit     => spinorbit_diffusion_equation      ! Defines the Usadel diffusion equation (spin-orbit terms)
    procedure                 :: interface_spinorbit_a   => spinorbit_interface_equation_a    ! Defines the left  boundary condition  (spin-orbit terms)
    procedure                 :: interface_spinorbit_b   => spinorbit_interface_equation_b    ! Defines the right boundary condition  (spin-orbit terms)

    ! These methods contain the equations that describe spin-active interfaces
    procedure                 :: interface_spinactive_a  => spinactive_interface_equation_a   ! Defines the left  boundary condition (spin-active terms) [TODO]
    procedure                 :: interface_spinactive_b  => spinactive_interface_equation_b   ! Defines the right boundary condition (spin-active terms) [TODO]

    ! These methods define miscellaneous utility functions
    procedure                 :: save                    => conductor_save                    ! Saves the state of the conductor to a different object
    procedure                 :: load                    => conductor_load                    ! Loads the state of the conductor from a different object
    procedure                 :: write_dos               => conductor_write_dos               ! Writes the density of states to a given output unit

    ! These methods are used by internal subroutines 
    final                     ::                            conductor_destruct                ! Type destructor
  end type

  ! Type constructors
  interface conductor
    module procedure conductor_construct
  end interface
contains

  !--------------------------------------------------------------------------------!
  !                IMPLEMENTATION OF CONSTRUCTORS AND DESTRUCTORS                  !
  !--------------------------------------------------------------------------------!

  pure function conductor_construct(energy, gap, thouless, scattering, points, spinorbit) result(this)
    ! Constructs a conductor object initialized to a superconducting state.
    type(conductor)                   :: this         ! Conductor object that will be constructed
    real(dp),    intent(in)           :: energy(:)    ! Discretized energy domain that will be used
    real(dp),    intent(in), optional :: thouless     ! Thouless energy       (default: see type declaration)
    real(dp),    intent(in), optional :: scattering   ! Imaginary energy term (default: see type declaration)
    complex(dp), intent(in), optional :: gap          ! Superconducting gap   (default: see definition below)
    integer,     intent(in), optional :: points       ! Number of positions   (default: see definition below)
    type(spin),  intent(in), optional :: spinorbit(:) ! Spin-orbit coupling   (default: zero coupling)
    integer                           :: n, m         ! Loop variables

    ! Optional argument: Thouless energy
    if (present(thouless)) then
      this%thouless = thouless
    end if

    ! Optional argument: imaginary energy
    if (present(scattering)) then
      this%scattering = scattering
    end if

    ! Optional argument: spin-orbit coupling
    if (allocated(this%spinorbit)) then
      deallocate(this%spinorbit)
    end if
    if (present(spinorbit)) then
      if (sum(spinorbit%norm()) > 1e-10) then
        allocate(this%spinorbit(size(spinorbit)))
        this%spinorbit = spinorbit
      end if
    end if
    
    ! Allocate memory (if necessary)
    if (.not. allocated(this%greenr)) then
      if (present(points)) then
        allocate(this%greenr(size(energy), points))
        allocate(this%energy(size(energy)))
        allocate(this%location(points))
      else
        allocate(this%greenr(size(energy), 150))
        allocate(this%energy(size(energy)))
        allocate(this%location(150))
      end if
    end if

    ! Initialize energy and position arrays
    this%energy   = energy
    this%location = [ ((real(n,kind=dp)/real(size(this%location)-1,kind=dp)), n=0,size(this%location)-1) ]

    ! Initialize the state
    if (present(gap)) then
      call this%init( gap )
    else
      call this%init( cmplx(1.0_dp,0.0_dp,kind=dp) )
    end if
  end function

  pure subroutine conductor_destruct(this)
    ! Define the type destructor.
    type(conductor), intent(inout) :: this

    ! Deallocate memory (if necessary)
    if(allocated(this%greenr)) then
      deallocate(this%greenr)
      deallocate(this%location)
      deallocate(this%energy)
    end if

    if (allocated(this%spinorbit)) then
      deallocate(this%spinorbit)
    end if
  end subroutine

  pure subroutine conductor_init(this, gap)
    ! Define the default initializer.
    class(conductor), intent(inout) :: this
    complex(dp),      intent(in   ) :: gap
    integer                         :: n, m

    do m = 1,size(this%location)
      do n = 1,size(this%energy)
        this%greenr(n,m) = green( cmplx(this%energy(n),this%scattering,kind=dp), gap )
      end do
    end do
  end subroutine

  !--------------------------------------------------------------------------------!
  !                     IMPLEMENTATION OF CONDUCTOR EQUATIONS                      !
  !--------------------------------------------------------------------------------!

  pure subroutine conductor_diffusion_equation(this, e, z, g, gt, dg, dgt, d2g, d2gt)
    ! Use the diffusion equation to calculate the second-derivatives of the Riccati parameters at energy e and point z.
    class(conductor), intent(in   ) :: this
    complex(dp),      intent(in   ) :: e
    real(dp),         intent(in   ) :: z
    type(spin),       intent(in   ) :: g, gt, dg, dgt
    type(spin),       intent(inout) :: d2g, d2gt
    type(spin)                      :: N, Nt

    ! Calculate the normalization matrices
    N   = spin_inv( pauli0 - g*gt )
    Nt  = spin_inv( pauli0 - gt*g )

    ! Calculate the second-derivatives of the Riccati parameters
    d2g  = (-2.0_dp,0.0_dp)*dg*Nt*gt*dg - (0.0_dp,2.0_dp)*e*g
    d2gt = (-2.0_dp,0.0_dp)*dgt*N*g*dgt - (0.0_dp,2.0_dp)*e*gt

    ! Calculate the contribution from a spin-orbit coupling
    if (allocated(this%spinorbit)) then
      call this%diffusion_spinorbit(g, gt, dg, dgt, d2g, d2gt)
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
        else
          ! Interface is a tunnel junction
          call this%interface_tunnel_a(a, g, gt, dg, dgt, r, rt)
        end if
      else
        ! Interface is a vacuum junction
        call this%interface_vacuum_a(g, gt, dg, dgt, r, rt)
      end if

      if (allocated(this%spinorbit)) then
        ! Interface has spin-orbit coupling
        call this%interface_spinorbit_a(g, gt, dg, dgt, r, rt)
      end if

      if (allocated(this%magnetization_a)) then
        ! Interface is spin-active
        call this%interface_spinactive_a(a, g, gt, dg, dgt, r, rt)
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
        else
          ! Interface is a tunnel junction
          call this%interface_tunnel_b(b, g, gt, dg, dgt, r, rt)
        end if
      else
        ! Interface is a vacuum junction
        call this%interface_vacuum_b(g, gt, dg, dgt, r, rt)
      end if

      if (allocated(this%spinorbit)) then
        ! Interface has a spin-orbit coupling
        call this%interface_spinorbit_b(g, gt, dg, dgt, r, rt)
      end if

      if (allocated(this%magnetization_b)) then
        ! Interface is spin-active
        call this%interface_spinactive_b(b, g, gt, dg, dgt, r, rt)
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

    ! Prepare variables associated with spin-active interfaces
    call spinactive_update_prehook(this)

    ! Modify the type string
    if (colors) then
      this%type_string = color_yellow // 'CONDUCTOR' // color_none
      if (allocated(this%spinorbit))       this%type_string = trim(this%type_string) // color_cyan   // ' [SOC]' // color_none
      if (allocated(this%magnetization_a)) this%type_string = trim(this%type_string) // color_purple // ' [SAL]' // color_none
      if (allocated(this%magnetization_b)) this%type_string = trim(this%type_string) // color_purple // ' [SAR]' // color_none
    else
      this%type_string = 'CONDUCTOR'
      if (allocated(this%spinorbit))       this%type_string = trim(this%type_string) // ' [SOC]'
      if (allocated(this%magnetization_a)) this%type_string = trim(this%type_string) // ' [SAL]'
      if (allocated(this%magnetization_b)) this%type_string = trim(this%type_string) // ' [SAR]'
    end if
  end subroutine

  impure subroutine conductor_update_posthook(this)
    ! Code to execute after running the update method of a class(conductor) object.
    class(conductor), intent(inout) :: this

    continue
  end subroutine

  !--------------------------------------------------------------------------------!
  !                    IMPLEMENTATION OF INPUT/OUTPUT METHODS                      !
  !--------------------------------------------------------------------------------!

  impure subroutine conductor_save(this, other)
    ! Saves the state of the conductor to a different object.
    class(conductor), intent(inout) :: this
    class(conductor), intent(inout) :: other

    other % greenr = this % greenr
  end subroutine

  impure subroutine conductor_load(this, other)
    ! Saves the state of the conductor to a different object.
    class(conductor), intent(inout) :: this
    class(conductor), intent(inout) :: other

    this % greenr = other % greenr
  end subroutine

  impure subroutine conductor_write_dos(this, unit, a, b)
    ! Writes the density of states as a function of position and energy to a given output unit.
    class(conductor),   intent(in) :: this      ! Material that the density of states will be calculated from
    integer,            intent(in) :: unit      ! Output unit that determines where the information will be written
    real(dp),           intent(in) :: a, b      ! Left and right end points of the material
    integer                        :: n, m      ! Temporary loop variables

    if (minval(this%energy) < -1e-16_dp) then
      ! If we have data for both positive and negative energies, simply write out the data
      do m=1,size(this%location)
        do n=1,size(this%energy)
          write(unit,*) a+(b-a)*this%location(m), this%energy(n), this%greenr(n,m)%get_dos()
        end do
        write(unit,*)
      end do
    else
      ! If we only have data for positive energies, assume that the negative region is symmetric
      do m=1,size(this%location)
        do n=size(this%energy),1,-1
          write(unit,*) a+(b-a)*this%location(m), -this%energy(n), this%greenr(n,m)%get_dos()
        end do
        do n=1,size(this%energy),+1
          write(unit,*) a+(b-a)*this%location(m), +this%energy(n), this%greenr(n,m)%get_dos()
        end do
        write(unit,*)
      end do
    end if
  end subroutine

  !--------------------------------------------------------------------------------!
  !                     IMPLEMENTATION OF SPIN-ORBIT COUPLING                      !
  !--------------------------------------------------------------------------------!

  ! TODO: These methods should be moved to a submodule when GFortran supports that.

  pure subroutine spinorbit_update_prehook(this)
    ! Updates the internal variables associated with spin-orbit coupling.
    class(conductor), intent(inout) :: this 

    if (allocated(this%spinorbit)) then
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
  end subroutine

  pure subroutine spinorbit_diffusion_equation(this, g, gt, dg, dgt, d2g, d2gt)
    ! Calculate the spin-orbit coupling terms in the diffusion equation, and update the second derivatives of the Riccati parameters.
    class(conductor), target, intent(in   ) :: this
    type(spin),               intent(in   ) :: g, gt, dg, dgt
    type(spin),               intent(inout) :: d2g, d2gt
    type(spin)                              :: N,  Nt

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
         + (2.0_dp,0.0_dp) * (Ax * g + g * Axt) * Nt * (Axt + gt * Ax * g)  &
         + (2.0_dp,0.0_dp) * (Ay * g + g * Ayt) * Nt * (Ayt + gt * Ay * g)  &
         + (2.0_dp,0.0_dp) * (Az * g + g * Azt) * Nt * (Azt + gt * Az * g)  &
         + (0.0_dp,2.0_dp) * (Az + g * Azt * gt) * N * dg                   &
         + (0.0_dp,2.0_dp) * dg * Nt * (gt * Az * g + Azt)

    d2gt = d2gt            + (A2t * gt - gt * A2)                           &
         + (2.0_dp,0.0_dp) * (Axt * gt + gt * Ax) * N * (Ax + g * Axt * gt) &
         + (2.0_dp,0.0_dp) * (Ayt * gt + gt * Ay) * N * (Ay + g * Ayt * gt) &
         + (2.0_dp,0.0_dp) * (Azt * gt + gt * Az) * N * (Az + g * Azt * gt) &
         - (0.0_dp,2.0_dp) * (Azt + gt * Az * g) * Nt * dgt                 &
         - (0.0_dp,2.0_dp) * dgt * N * (g * Azt * gt + Az)

    end associate
  end subroutine

  pure subroutine spinorbit_interface_equation_a(this, g1, gt1, dg1, dgt1, r1, rt1)
    ! Calculate the spin-orbit coupling terms in the left boundary condition, and update the residuals.
    class(conductor), target, intent(in   ) :: this
    type(spin),               intent(in   ) :: g1, gt1, dg1, dgt1
    type(spin),               intent(inout) :: r1, rt1

    ! Rename the spin-orbit coupling matrices
    associate(Az  => this % Az,&
              Azt => this % Azt)

    ! Update the residuals
    r1  = r1  - (0.0_dp,1.0_dp) * (Az  * g1  + g1  * Azt)
    rt1 = rt1 + (0.0_dp,1.0_dp) * (Azt * gt1 + gt1 * Az )

    end associate
  end subroutine

  pure subroutine spinorbit_interface_equation_b(this, g2, gt2, dg2, dgt2, r2, rt2)
    ! Calculate the spin-orbit coupling terms in the right boundary condition, and update the residuals.
    class(conductor), target, intent(in   ) :: this
    type(spin),               intent(in   ) :: g2, gt2, dg2, dgt2
    type(spin),               intent(inout) :: r2, rt2

    ! Rename the spin-orbit coupling matrices
    associate(Az   => this % Az,&
              Azt  => this % Azt)

    ! Update the residuals
    r2  = r2  - (0.0_dp,1.0_dp) * (Az  * g2  + g2  * Azt)
    rt2 = rt2 + (0.0_dp,1.0_dp) * (Azt * gt2 + gt2 * Az )  

    end associate
  end subroutine

  !--------------------------------------------------------------------------------!
  !                   IMPLEMENTATION OF SPIN-ACTIVE INTERFACES                     !
  !--------------------------------------------------------------------------------!

  ! TODO: These methods should be moved to a submodule when GFortran supports that.

  pure subroutine spinactive_update_prehook(this)
    ! Updates the internal variables associated with spin-active interfaces.
    class(conductor), intent(inout) :: this 

    if (allocated(this % magnetization_a)) then
      ! Rename the relevant variables
      associate(G0 => this % conductance_a,   &
                G1 => this % G1_a,            &
                GM => this % GM_a,            &
                GP => this % GP_a,            &
                h  => this % magnetization_a, & 
                P  => this % polarization_a,  &
                Q  => this % phaseshift_a,    &
                M  => this % M_a)

      ! Calculate the elements of the magnetization matrix diag(m·σ,m·σ*)
      M(1:2,1:2) = h(1) * pauli1 + h(2) * pauli2 + h(3) * pauli3
      M(3:4,3:4) = h(1) * pauli1 - h(2) * pauli2 + h(3) * pauli3

      ! Calculate the conductances associated with the spin-active interface
      GM = G0 * P/(1 + sqrt(1-P**2))
      GP = G0 * (-2*i*Q)/(1 + sqrt(1-P**2))
      G1 = G0 * (1 - sqrt(1-P**2))/(1 + sqrt(1-P**2))

      end associate
    end if

    if (allocated(this % magnetization_b)) then
      ! Rename the relevant variables
      associate(G0 => this % conductance_b,   &
                G1 => this % G1_b,            &
                GM => this % GM_b,            &
                GP => this % GP_b,            &
                h  => this % magnetization_b, & 
                P  => this % polarization_b,  &
                Q  => this % phaseshift_b,    &
                M  => this % M_b)

      ! Calculate the elements of the magnetization matrix diag(m·σ,m·σ*)
      M(1:2,1:2) = h(1) * pauli1 + h(2) * pauli2 + h(3) * pauli3
      M(3:4,3:4) = h(1) * pauli1 - h(2) * pauli2 + h(3) * pauli3

      ! Calculate the conductances associated with the spin-active interface
      GM = G0 * P/(1 + sqrt(1-P**2))
      GP = G0 * (-2*i*Q)/(1 + sqrt(1-P**2))
      G1 = G0 * (1 - sqrt(1-P**2))/(1 + sqrt(1-P**2))

      end associate
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
    associate(G1 => this % G1_a,&
              GM => this % GM_a,&
              GP => this % GP_a,&
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

    ! Calculate the 8×8 Green's function in the left material
    L(1:2,1:2) = (+1.0_dp) * N0  * (pauli0 + g0*gt0)
    L(1:2,3:4) = (+2.0_dp) * N0  * g0
    L(3:4,1:2) = (-2.0_dp) * Nt0 * gt0
    L(3:4,3:4) = (-1.0_dp) * Nt0 * (pauli0 + gt0*g0)

    ! Calculate the 8×8 Green's function in the right material
    R(1:2,1:2) = (+1.0_dp) * N1  * (pauli0 + g1*gt1)
    R(1:2,3:4) = (+2.0_dp) * N1  * g1
    R(3:4,1:2) = (-2.0_dp) * Nt1 * gt1
    R(3:4,3:4) = (-1.0_dp) * Nt1 * (pauli0 + gt1*g1)

    ! Calculate the spin-active terms in the interface current
    I = G1 * (matmul(R,matmul(M,matmul(L,M)))  &
             -matmul(M,matmul(L,matmul(M,R)))) &
      + GM * (matmul(R,matmul(L,M)+matmul(M,L))        &
             -matmul(matmul(L,M)+matmul(M,L),R))       &
      + GP * (matmul(R,M) - matmul(M,R))

    ! Calculate the deviation from the boundary condition
    r1  = r1  - (pauli0 - g1*gt1) * (I(1:2,3:4) - I(1:2,1:2)*g1)
    rt1 = rt1 - (pauli0 - gt1*g1) * (I(3:4,1:2) - I(3:4,3:4)*gt1)

    end associate
    end associate
  end subroutine

  pure subroutine spinactive_interface_equation_b(this, b, g2, gt2, dg2, dgt2, r2, rt2)
    ! Calculate the spin-active terms in the right boundary condition, and update the residuals.
    class(conductor), target, intent(in   ) :: this
    type(green),              intent(in   ) :: b
    type(spin),               intent(in   ) :: g2, gt2, dg2, dgt2
    type(spin),               intent(inout) :: r2, rt2

    ! TODO
    continue
  end subroutine
end module
