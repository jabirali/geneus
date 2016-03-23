! This module defines the data type 'halfmetal', which models the physical state of a strong or halfmetallic ferromagnet.
! @TODO: Finish this module.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2016-03-08
! Updated: 2016-03-08

module mod_halfmetal
  use mod_stdio
  use mod_math
  use mod_spin
  use mod_green
  use mod_material
  implicit none
  private

  ! Type declarations
  type, public, extends(material) :: halfmetal
    ! These parameters control the physical characteristics of the material
    real(wp)                  :: polarization            =  0.99_wp                        ! Spin-polarization of the ferromagnet
  contains
    ! These methods are required by the class(material) abstract interface
    procedure                 :: init                    => halfmetal_init                 ! Initializes the propagators
    procedure                 :: diffusion_equation      => halfmetal_diffusion_equation   ! Defines the Usadel diffusion equation
    procedure                 :: interface_equation_a    => halfmetal_interface_equation_a ! Boundary condition at the left  interface
    procedure                 :: interface_equation_b    => halfmetal_interface_equation_b ! Boundary condition at the right interface
    procedure                 :: update_prehook          => halfmetal_update_prehook       ! Code to execute before calculating the propagators
    procedure                 :: update_posthook         => halfmetal_update_posthook      ! Code to execute after  calculating the propagators

    ! These methods define miscellaneous utility functions
    procedure                 :: conf                    => halfmetal_conf                 ! Configures material parameters
  end type

  ! Type constructors
  interface halfmetal
    module procedure halfmetal_construct
  end interface
contains

  !--------------------------------------------------------------------------------!
  !                        IMPLEMENTATION OF CONSTRUCTORS                          !
  !--------------------------------------------------------------------------------!

  pure function halfmetal_construct(cutoff) result(this)
    ! Constructs a halfmetal object.
    type(halfmetal)                :: this   ! Halfmetal object that will be constructed
    real(wp), intent(in), optional :: cutoff ! Debye cutoff for the energy domain

    ! Initialize locations
    allocate(this%location(150))
    call linspace(this%location, 0.0_wp, 1.0_wp)

    ! Initialize energies
    allocate(this%energy(600))
    call linspace(this%energy(   :400), 1e-6_wp, 1.50_wp)
    call linspace(this%energy(400:500), 1.50_wp, 4.50_wp)
    call linspace(this%energy(500:   ), 4.50_wp, 30.0_wp)

    ! Initialize propagators
    allocate(this%propagator(size(this%energy),size(this%location)))
    call this%init( cx(1.0_wp,0.0_wp) )
  end function

  pure subroutine halfmetal_init(this, gap)
    ! Define the default initializer.
    class(halfmetal), intent(inout) :: this
    complex(wp),      intent(in   ) :: gap
    integer                         :: n, m

    do m = 1,size(this%location)
      do n = 1,size(this%energy)
        this%propagator(n,m) = green0
      end do
    end do
  end subroutine

  !--------------------------------------------------------------------------------!
  !                     IMPLEMENTATION OF HALFMETAL EQUATIONS                      !
  !--------------------------------------------------------------------------------!

  pure subroutine halfmetal_diffusion_equation(this, e, z, g, gt, dg, dgt, d2g, d2gt)
    ! Use the diffusion equation to calculate the second-derivatives of the Riccati parameters at energy e and point z.
    class(halfmetal), intent(in   ) :: this
    complex(wp),      intent(in   ) :: e
    real(wp),         intent(in   ) :: z
    type(spin),       intent(inout) :: d2g, d2gt
    type(spin),       intent(in   ) :: g, gt, dg, dgt
    type(spin)                      :: h, ht, dh, dht
    type(spin)                      :: N, Nt
    type(spin)                      :: P

    ! Ensure that the Riccati parameters are diagonal
    h   % matrix(1,2) = 0
    h   % matrix(2,1) = 0
    ht  % matrix(1,2) = 0
    ht  % matrix(2,1) = 0
    dh  % matrix(1,2) = 0
    dh  % matrix(2,1) = 0
    dht % matrix(1,2) = 0
    dht % matrix(2,1) = 0

    ! Calculate the polarization matrix
    P   % matrix(1,1) = 2/(1 + eps + this%polarization)
    P   % matrix(1,2) = 0
    P   % matrix(2,1) = 0
    P   % matrix(2,2) = 2/(1 + eps - this%polarization)

    ! Calculate the normalization matrices
    N   = spin_inv( pauli0 - h*ht )
    Nt  = spin_inv( pauli0 - ht*h )

    ! Calculate the second-derivatives of the Riccati parameters
    d2g  = (-2.0_wp,0.0_wp)*dh*Nt*ht*dh - (0.0_wp,2.0_wp)*e*P*h
    d2gt = (-2.0_wp,0.0_wp)*dht*N*h*dht - (0.0_wp,2.0_wp)*e*P*ht
  end subroutine

  pure subroutine halfmetal_interface_equation_a(this, a, g, gt, dg, dgt, r, rt)
      ! Calculate residuals from the boundary conditions at the left interface.
      class(halfmetal),          intent(in   ) :: this
      type(green),               intent(in   ) :: a
      type(spin),                intent(in   ) :: g, gt, dg, dgt
      type(spin),                intent(inout) :: r, rt

    ! Calculate the deviation between the materials
    associate(g0  => a%g, gt0 => a%gt)
      r  = g  - g0
      rt = gt - gt0
    end associate
  end subroutine

  pure subroutine halfmetal_interface_equation_b(this, b, g, gt, dg, dgt, r, rt)
    ! Calculate residuals from the boundary conditions at the right interface.
    class(halfmetal),          intent(in   ) :: this
    type(green),               intent(in   ) :: b
    type(spin),                intent(in   ) :: g, gt, dg, dgt
    type(spin),                intent(inout) :: r, rt

    ! Calculate the deviation between the materials
    associate(g3 => b%g, gt3 => b%gt)
      r  = g  - g3
      rt = gt - gt3
    end associate
  end subroutine

  impure subroutine halfmetal_update_prehook(this)
    ! Code to execute before running the update method of a class(halfmetal) object.
    class(halfmetal), intent(inout) :: this

    ! Modify the type string
    this%type_string = color_purple // 'HALFMETAL' // color_none
  end subroutine

  impure subroutine halfmetal_update_posthook(this)
    ! Code to execute after running the update method of a class(halfmetal) object.
    class(halfmetal), intent(inout) :: this

    ! Calculate the density of states
    call this%update_density

    ! Calculate the charge and spin currents
    call this%update_current
  end subroutine

  !--------------------------------------------------------------------------------!
  !                      IMPLEMENTATION OF UTILITY METHODS                         !
  !--------------------------------------------------------------------------------!

  impure subroutine halfmetal_conf(this, key, val)
    !! Configure a material property based on a key-value pair.
    class(halfmetal), intent(inout) :: this
    character(*),     intent(in   ) :: key
    character(*),     intent(in   ) :: val

    select case(key)
      case ('polarization')
        read(val,*) this % polarization
      case default
        call material_conf(this, key, val)
    end select
  end subroutine
end module
