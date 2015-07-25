! This module defines the data type 'conductor', which models the physical state of a conductor for a discretized range
! of positions and energies.  It has two main applications: (i) it can be used as a base type for more exotic materials,
! such as superconductors and ferromagnets; (ii) it can be used in conjunction with such materials in mesoscopic hybrid
! structures. The module also defines helper functions like 'connect' which operate on general class(conductor) objects.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-11
! Updated: 2015-07-25

module mod_conductor
  use mod_system
  use mod_spin
  use mod_green
  implicit none

  ! Type declarations
  type conductor
    ! These parameters control the physical characteristics of the material (should be modified by the user)
    real(dp)                  :: thouless      =  1.00_dp                           ! Thouless energy of the material (ratio of the diffusion constant to the squared material length)
    real(dp)                  :: scattering    =  0.01_dp                           ! Imaginary energy term (this models inelastic scattering processes and stabilizes the BVP solver)
    real(dp)                  :: conductance_a =  0.00_dp                           ! Tunneling conductance of the left interface  (relative to the bulk conductance of this material)
    real(dp)                  :: conductance_b =  0.00_dp                           ! Tunneling conductance of the right interface (relative to the bulk conductance of this material)
    class(conductor), pointer :: material_a    => null()                            ! Material connected to this one at the left  interface (default: null pointer, meaning vacuum)
    class(conductor), pointer :: material_b    => null()                            ! Material connected to this one at the right interface (default: null pointer, meaning vacuum)
    type(spin),   allocatable :: spinorbit(:)                                       ! This is an SU(2) vector field that describes spin-orbit coupling (default: no spin-orbit coupling)
  
    ! These parameters control the boundary value problem solver (can be modified by the user)
    integer                   :: scaling       =  64                                ! Maximal allowed scaling of the mesh resolution (range: 2^N, N>1)
    integer                   :: order         =  4                                 ! Order of the Runge—Kutta method used by the solver (range: 2, 4, 6)
    integer                   :: control       =  2                                 ! Error control method (1: defect, 2: global error, 3: 1 then 2, 4: 1 and 2)
    integer                   :: information   =  0                                 ! Information amount that should be written to standard out (range: [-1,2])
    real(dp)                  :: tolerance     =  1e-4_dp                           ! Error tolerance (determines the maximum allowed defect or error)

    ! These variables store the physical state of the material (should not be modified by the user)
    type(green), allocatable  :: state(:,:)                                         ! Physical state as a function of energy and position
    real(dp),    allocatable  :: energy(:)                                          ! Discretized  energy  domain that will be considered
    real(dp),    allocatable  :: location(:)                                        ! Discretized position domain that will be considered

    ! These variables are used by internal subroutines (should not be accessed by the user)
    complex(dp)               :: erg                                                ! Temporary storage for the current energy
    type(spin)                :: Ax,  Ay,  Az,  A2                                  ! Spin-orbit coupling matrices (Cartesian components and square)
    type(spin)                :: Axt, Ayt, Azt, A2t                                 ! Spin-orbit coupling matrices (tilde-conjugated versions)
    type(green), pointer      :: state_a       => null()                            ! Pointer to the left  interface state for the current energy
    type(green), pointer      :: state_b       => null()                            ! Pointer to the right interface state for the current energy

    ! Miscellaneous variables
    character(len=64)         :: type_string   =  ''                                ! This variable should be modified when class(conductor) subtypes are initialized
  contains
    ! These methods control the simulation process (can be invoked by the user)
    procedure                 :: initialize         => conductor_initialize         ! Resets  the internal state of the material
    procedure                 :: update             => conductor_update             ! Updates the internal state of the material

    ! These methods contain the equations that describe the material (should not be invoked by the user)
    procedure                 :: usadel_equation    => conductor_usadel_equation    ! Differential equation that describes the conductor
    procedure                 :: interface_vacuum_a => conductor_interface_vacuum_a ! Defines the left  boundary condition for a vacuum interface
    procedure                 :: interface_vacuum_b => conductor_interface_vacuum_b ! Defines the right boundary condition for a vacuum interface
    procedure                 :: interface_tunnel_a => conductor_interface_tunnel_a ! Defines the left  boundary condition for a tunnel interface
    procedure                 :: interface_tunnel_b => conductor_interface_tunnel_b ! Defines the right boundary condition for a tunnel interface
    procedure                 :: update_fields      => conductor_update_fields      ! Updates the physical fields based on stored Green's functions

    ! These methods are used to output physical results (should be invoked by the user)
    procedure                 :: write_dos          => conductor_write_dos          ! Writes the density of states to a given output unit

    ! These methods are used by internal subroutines (should not be invoked by the user)
    final                     ::                       conductor_destruct           ! Type destructor
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
    if (present(spinorbit)) then
      allocate(this%spinorbit(size(spinorbit)))
      this%spinorbit = spinorbit
    end if
    
    ! Allocate memory (if necessary)
    if (.not. allocated(this%state)) then
      if (present(points)) then
        allocate(this%state(size(energy), points))
        allocate(this%location(points))
        allocate(this%energy(size(energy)))
      else
        allocate(this%state(size(energy), 150))
        allocate(this%location(150))
        allocate(this%energy(size(energy)))
      end if
    end if

    ! Initialize energy and position arrays
    this%energy   = energy
    this%location = [ ((real(n,kind=dp)/real(size(this%location)-1,kind=dp)), n=0,size(this%location)-1) ]

    ! Initialize the state
    if (present(gap)) then
      call this%initialize( gap )
    else
      call this%initialize( cmplx(1.0_dp,0.0_dp,kind=dp) )
    end if

    ! Modify the type string
    this%type_string = color_yellow // 'CONDUCTOR' // color_none
    if (allocated(this%spinorbit)) then
      this%type_string = trim(this%type_string) // color_cyan // ' [SOC] ' // color_none
    end if
  end function

  pure subroutine conductor_destruct(this)
    ! Define the type destructor.
    type(conductor), intent(inout) :: this

    ! Deallocate memory (if necessary)
    if(allocated(this%state)) then
      deallocate(this%state)
      deallocate(this%location)
      deallocate(this%energy)
    end if

    if (allocated(this%spinorbit)) then
      deallocate(this%spinorbit)
    end if
  end subroutine

  pure subroutine conductor_initialize(this, gap)
    ! Define the default initializer.
    class(conductor), intent(inout) :: this
    complex(dp),      intent(in   ) :: gap
    integer                         :: n, m

    do m = 1,size(this%location)
      do n = 1,size(this%energy)
        this%state(n,m) = green( cmplx(this%energy(n),this%scattering,kind=dp), gap )
      end do
    end do
  end subroutine

  !--------------------------------------------------------------------------------!
  !                     IMPLEMENTATION OF BVP SOLVER METHODS                       !
  !--------------------------------------------------------------------------------!

  subroutine conductor_update(this)
    ! This subroutine updates the current estimate for the state of the system.
    use bvp_m

    class(conductor), intent(inout) :: this                       ! Object that will be updated
    real(dp)                        :: u(32,size(this%location))  ! Real state vector required by the BVP solver
    type(bvp_sol)                   :: sol                        ! Object with information about the BVP solution
    integer                         :: n, m                       ! Internal loop variables

    ! Prepare variables associated with spin-orbit coupling
    if (allocated(this%spinorbit)) then
      call spinorbit_update(this)
    end if

    ! Status information
    if (this%information >= 0) then
      write(stdout,'(a)') color_white // ' :: ' // color_none // trim(this%type_string) // '                                     '
    end if

    do n=1,size(this%energy)
      ! Status information
      if (this%information >= 0) then
        write(stdout,'(4x,a,1x,i4,1x,a,1x,i4,1x,a,f0.5,a1)',advance='no') &
          '[',n,'/',size(this%energy),']  ϵ = ',this%energy(n), achar(13)
        flush(stdout)
      end if

      ! Convert all states at this energy level to real-valued state vectors
      do m=1,size(this%location)
        u(:,m) = this%state(n,m)
      end do

      ! Copy the complex energy and boundary conditions to internal work variables
      this%erg = cmplx(this%energy(n)/this%thouless, this%scattering/this%thouless, kind=dp)
      if (associated(this%material_a)) then
        this%state_a => this%material_a%state(n,ubound(this%material_a%state,2))
      end if
      if (associated(this%material_b)) then
        this%state_b => this%material_b%state(n,lbound(this%material_b%state,2))
      end if

      ! Initialize the BVP solver
      sol = bvp_init(32, 16, this%location, u, max_num_subintervals=(size(this%location)*this%scaling))

      ! Solve the differential equation
      sol = bvp_solver(sol, ode, bc, method=this%order, error_control=this%control, tol=this%tolerance, trace=this%information)

      ! Use the results to update the state
      call bvp_eval(sol, this%location, u)
      do m=1,size(this%location)
        this%state(n,m) = u(:,m)
      end do
    end do

    ! Clean up
    call bvp_terminate(sol)

    ! Update other internal variables if necessary
    call this%update_fields   
  contains
    subroutine ode(z, u, f)
      ! Definition of the differential equation u'=f(z,u).
      real(dp), intent(in)  :: z
      real(dp), intent(in)  :: u(32)
      real(dp), intent(out) :: f(32)

      type(spin)            :: g, gt, dg, dgt, d2g, d2gt

      ! Extract the Riccati parameters
      g   = u( 1: 8)
      gt  = u( 9:16)
      dg  = u(17:24)
      dgt = u(25:32)

      ! Calculate the second derivatives of the Riccati parameters
      call this%usadel_equation(z, g, gt, dg, dgt, d2g, d2gt)

      ! Calculate the contribution from a spin-orbit coupling
      if (allocated(this%spinorbit)) then
        call spinorbit_usadel_equation(this, g, gt, dg, dgt, d2g, d2gt)
      end if
            
      ! Pack the results into a state vector
      f( 1: 8) = dg
      f( 9:16) = dgt
      f(17:24) = d2g
      f(25:32) = d2gt
    end subroutine

    subroutine bc(ua, ub, bca, bcb)
      ! Definition of the boundary conditions bca=g(ua) and bcb=g(ub).
      real(dp), intent(in)  :: ua(32)
      real(dp), intent(in)  :: ub(32)
      real(dp), intent(out) :: bca(16)
      real(dp), intent(out) :: bcb(16)

      type(spin)            :: g1, gt1, dg1, dgt1, r1, rt1
      type(spin)            :: g2, gt2, dg2, dgt2, r2, rt2

      ! State at the left end of the material
      g1   = ua( 1: 8)
      gt1  = ua( 9:16)
      dg1  = ua(17:24)
      dgt1 = ua(25:32)

      ! State at the right end of the material
      g2   = ub( 1: 8)
      gt2  = ub( 9:16)
      dg2  = ub(17:24)
      dgt2 = ub(25:32)

      ! Calculate residuals from the boundary conditions
      if (associated(this%material_a)) then
        ! Left interface is a tunnel junction
        call this%interface_tunnel_a(g1, gt1, dg1, dgt1, r1, rt1)
      else
        ! Left interface is a vacuum junction
        call this%interface_vacuum_a(g1, gt1, dg1, dgt1, r1, rt1)
      end if

      if (associated(this%material_b)) then
        ! Right interface is a tunnel junction
        call this%interface_tunnel_b(g2, gt2, dg2, dgt2, r2, rt2)
      else
        ! Right interface is a vacuum junction
        call this%interface_vacuum_b(g2, gt2, dg2, dgt2, r2, rt2)
      end if

      if (allocated(this%spinorbit)) then
        ! Boundary conditions must include spin-orbit coupling
        call spinorbit_interface_a(this, g1, gt1, dg1, dgt1, r1, rt1)
        call spinorbit_interface_b(this, g2, gt2, dg2, dgt2, r2, rt2)
      end if

      ! Pack the results into state vectors
      bca(1: 8) = r1
      bca(9:16) = rt1
      bcb(1: 8) = r2
      bcb(9:16) = rt2
    end subroutine
  end subroutine

  !--------------------------------------------------------------------------------!
  !                     IMPLEMENTATION OF CONDUCTOR EQUATIONS                      !
  !--------------------------------------------------------------------------------!

  subroutine conductor_usadel_equation(this, z, g, gt, dg, dgt, d2g, d2gt)
    ! Use the Usadel equation to calculate the second derivatives of the Riccati parameters at point z.
    class(conductor), intent(in   ) :: this
    real(dp),         intent(in   ) :: z
    type(spin),       intent(in   ) :: g, gt, dg, dgt
    type(spin),       intent(inout) :: d2g, d2gt
    type(spin)                      :: N, Nt

    ! Calculate the normalization matrices
    N   = spin_inv( pauli0 - g*gt )
    Nt  = spin_inv( pauli0 - gt*g )

    ! Calculate the second derivatives of the Riccati parameters
    d2g  = (-2.0_dp,0.0_dp)*dg*Nt*gt*dg - (0.0_dp,2.0_dp)*this%erg*g
    d2gt = (-2.0_dp,0.0_dp)*dgt*N*g*dgt - (0.0_dp,2.0_dp)*this%erg*gt
  end subroutine

  subroutine conductor_interface_vacuum_a(this, g1, gt1, dg1, dgt1, r1, rt1)
    ! Defines a vacuum boundary condition for the left interface.
    class(conductor), intent(in   ) :: this
    type(spin),       intent(in   ) :: g1, gt1, dg1, dgt1
    type(spin),       intent(inout) :: r1, rt1

    r1  = dg1
    rt1 = dgt1
  end subroutine

  subroutine conductor_interface_vacuum_b(this, g2, gt2, dg2, dgt2, r2, rt2)
    ! Defines a vacuum boundary condition for the right interface.
    class(conductor), intent(in   ) :: this
    type(spin),       intent(in   ) :: g2, gt2, dg2, dgt2
    type(spin),       intent(inout) :: r2, rt2

    r2  = dg2
    rt2 = dgt2
  end subroutine

  subroutine conductor_interface_tunnel_a(this, g1, gt1, dg1, dgt1, r1, rt1)
    ! Defines a tunneling boundary condition for the left interface.
    class(conductor), intent(in   ) :: this
    type(spin),       intent(inout) :: r1, rt1
    type(spin),       intent(in   ) :: g1, gt1, dg1, dgt1
    type(spin),       pointer       :: g0, gt0, dg0, dgt0
    type(spin)                      :: N0, Nt0

    ! Rename the state in the material to the left
    g0   => this%state_a%g
    gt0  => this%state_a%gt
    dg0  => this%state_a%dg
    dgt0 => this%state_a%dgt

    ! Calculate the normalization matrices
    N0  = spin_inv( pauli0 - g0*gt0 )
    Nt0 = spin_inv( pauli0 - gt0*g0 )

    ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
    r1  = dg1  - this%conductance_a*( pauli0 - g1*gt0 )*N0*(  g1  - g0  )
    rt1 = dgt1 - this%conductance_a*( pauli0 - gt1*g0 )*Nt0*( gt1 - gt0 )
  end subroutine

  subroutine conductor_interface_tunnel_b(this, g2, gt2, dg2, dgt2, r2, rt2)
    ! Defines a tunneling boundary condition for the right interface.
    class(conductor), intent(in   ) :: this
    type(spin),       intent(inout) :: r2, rt2
    type(spin),       intent(in   ) :: g2, gt2, dg2, dgt2
    type(spin),       pointer       :: g3, gt3, dg3, dgt3
    type(spin)                      :: N3, Nt3

    ! Rename the state in the material to the right
    g3   => this%state_b%g
    gt3  => this%state_b%gt
    dg3  => this%state_b%dg
    dgt3 => this%state_b%dgt
  
    ! Calculate the normalization matrices
    N3  = spin_inv( pauli0 - g3*gt3 )
    Nt3 = spin_inv( pauli0 - gt3*g3 )

    ! Calculate the deviation from the Kuprianov--Lukichev boundary condition
    r2  = dg2  - this%conductance_b*( pauli0 - g2*gt3 )*N3*(  g3  - g2  )
    rt2 = dgt2 - this%conductance_b*( pauli0 - gt2*g3 )*Nt3*( gt3 - gt2 )
  end subroutine

  subroutine conductor_update_fields(this)
    ! Define an interface used by class(conductor) objects to self-consistently update physical fields.
    class(conductor), intent(inout) :: this

    continue
  end subroutine

  !--------------------------------------------------------------------------------!
  !                    IMPLEMENTATION OF INPUT/OUTPUT METHODS                      !
  !--------------------------------------------------------------------------------!

  subroutine conductor_write_dos(this, unit, a, b)
    ! Writes the density of states as a function of position and energy to a given output unit.
    class(conductor),   intent(in) :: this      ! Material that the density of states will be calculated from
    integer,            intent(in) :: unit      ! Output unit that determines where the information will be written
    real(dp),           intent(in) :: a, b      ! Left and right end points of the material
    integer                        :: n, m      ! Temporary loop variables

    if (minval(this%energy) < 0.0_dp) then
      ! If we have data for both positive and negative energies, simply write out the data
      do m=1,size(this%location)
        do n=1,size(this%energy)
          write(unit,*) a+(b-a)*this%location(m), this%energy(n), this%state(n,m)%get_dos()
        end do
        write(unit,*)
      end do
    else
      ! If we only have data for positive energies, assume that the negative region is symmetric
      do m=1,size(this%location)
        do n=size(this%energy),1,-1
          write(unit,*) a+(b-a)*this%location(m), -this%energy(n), this%state(n,m)%get_dos()
        end do
        do n=1,size(this%energy),+1
          write(unit,*) a+(b-a)*this%location(m), +this%energy(n), this%state(n,m)%get_dos()
        end do
        write(unit,*)
      end do
    end if
  end subroutine

  !--------------------------------------------------------------------------------!
  !                 IMPLEMENTATION OF SPIN-ORBIT COUPLING METHODS                  !
  !--------------------------------------------------------------------------------!

  subroutine spinorbit_update(this)
    ! Updates the internal variables associated with a spin-orbit coupling.
    class(conductor), intent(inout) :: this 

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
  end subroutine

  subroutine spinorbit_usadel_equation(this, g, gt, dg, dgt, d2g, d2gt)
    ! Calculate the spin-orbit coupling terms in the Usadel equation, and update the second derivatives of the Riccati parameters.
    type(conductor),  target, intent(in   ) :: this
    type(spin),               intent(in   ) :: g, gt, dg, dgt
    type(spin),               intent(inout) :: d2g, d2gt

    type(spin)                              :: N,  Nt
    type(spin),       pointer               :: Ax, Axt
    type(spin),       pointer               :: Ay, Ayt
    type(spin),       pointer               :: Az, Azt
    type(spin),       pointer               :: A2, A2t

    ! Calculate the normalization matrices
    N   = spin_inv( pauli0 - g*gt )
    Nt  = spin_inv( pauli0 - gt*g )

    ! Rename the spin-orbit coupling matrices
    Ax => this % Ax;  Axt => this % Axt;
    Ay => this % Ay;  Ayt => this % Ayt;
    Az => this % Az;  Azt => this % Azt;
    A2 => this % A2;  A2t => this % A2t;

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
  end subroutine

  subroutine spinorbit_interface_a(this, g1, gt1, dg1, dgt1, r1, rt1)
    ! Calculate the spin-orbit coupling terms in the left boundary condition, and update the residuals.
    type(conductor), target, intent(in   ) :: this
    type(spin),              intent(in   ) :: g1, gt1, dg1, dgt1
    type(spin),              intent(inout) :: r1, rt1
    type(spin),      pointer               :: Az, Azt

    ! Rename the spin-orbit coupling matrices
    Az   => this % Az
    Azt  => this % Azt

    ! Update the residuals
    r1  = r1  - (0.0_dp,1.0_dp) * (Az  * g1  + g1  * Azt)
    rt1 = rt1 + (0.0_dp,1.0_dp) * (Azt * gt1 + gt1 * Az )
  end subroutine

  subroutine spinorbit_interface_b(this, g2, gt2, dg2, dgt2, r2, rt2)
    ! Calculate the spin-orbit coupling terms in the right boundary condition, and update the residuals.
    type(conductor), target, intent(in   ) :: this
    type(spin),              intent(in   ) :: g2, gt2, dg2, dgt2
    type(spin),              intent(inout) :: r2, rt2
    type(spin),      pointer               :: Az, Azt

    ! Rename the spin-orbit coupling matrices
    Az   => this % Az
    Azt  => this % Azt

    ! Update the residuals
    r2  = r2  - (0.0_dp,1.0_dp) * (Az  * g2  + g2  * Azt)
    rt2 = rt2 + (0.0_dp,1.0_dp) * (Azt * gt2 + gt2 * Az )  
  end subroutine
end module
