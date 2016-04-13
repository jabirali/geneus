!> Author:   Jabir Ali Ouassou
!> Date:     2016-03-08
!> Category: Materials
!>
!> This module defines the data type 'halfmetal', which models the physical state of a strong or halfmetallic ferromagnet.
!> The type is a member of class(conductor), and inherits the internal structure and generic methods defined there.
!>
!> @TODO 
!>   Add an update_posthook to rescale density(:) and current(:) with the dependence of the diffusion constant matrix
!>   on the polarization. Remember to check how the polarization dependence varies with the number of dimensions.
!>
!> @TODO
!>   Check if a non-linear dependence of the polarization matrix on the polarization is more sensible?

module halfmetal_m
  use stdio_m
  use math_m
  use spin_m
  use green_m
  use material_m
  use conductor_m
  implicit none
  private

  ! Type declarations
  type, public, extends(conductor) :: halfmetal
    real(wp)            :: polarization =  0.00_wp                                ! Spin-polarization of the ferromagnet
    type(spin), private :: P                                                      ! Polarization matrix
  contains
    procedure           :: init                 => halfmetal_init                 ! Initializes the propagators
    procedure           :: conf                 => halfmetal_conf                 ! Configures the material parameters
    procedure           :: diffusion_equation   => halfmetal_diffusion_equation   ! Defines the Usadel diffusion equation
    procedure           :: interface_equation_a => halfmetal_interface_equation_a ! Boundary condition at the left  interface
    procedure           :: interface_equation_b => halfmetal_interface_equation_b ! Boundary condition at the right interface
    procedure           :: update_prehook       => halfmetal_update_prehook       ! Code to execute before calculating the propagators
    procedure           :: update_posthook      => halfmetal_update_posthook      ! Code to execute after  calculating the propagators
  end type

  ! Type constructors
  interface halfmetal
    module procedure halfmetal_construct
  end interface
contains

  !--------------------------------------------------------------------------------!
  !                        IMPLEMENTATION OF CONSTRUCTORS                          !
  !--------------------------------------------------------------------------------!

  function halfmetal_construct() result(this)
    ! Constructs a halfmetal object.
    type(halfmetal) :: this

    ! Call the superclass constructor
    this%conductor = conductor()

    ! Initialize the propagators
    call this%init(cx(0.0_wp))
  end function

  pure subroutine halfmetal_init(this, gap)
    ! Initializes the propagators to a non-superconducting state.
    class(halfmetal), intent(inout) :: this
    complex(wp),      intent(in   ) :: gap
    integer                         :: n, m

    do m = 1,size(this%location)
      do n = 1,size(this%energy)
        this%propagator(n,m) = green()
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

    ! Ensure that the Riccati parameters are diagonal
    h   = g   % matrix * pauli0 % matrix
    ht  = gt  % matrix * pauli0 % matrix
    dh  = dg  % matrix * pauli0 % matrix
    dht = dgt % matrix * pauli0 % matrix

    ! Calculate the normalization matrices
    N   = spin_inv( pauli0 - h*ht )
    Nt  = spin_inv( pauli0 - ht*h )

    ! Calculate the second-derivatives of the Riccati parameters
    associate(P => this % P)
      d2g  = (-2.0_wp,0.0_wp)*dh*Nt*ht*dh - (0.0_wp,2.0_wp)*e*P*h
      d2gt = (-2.0_wp,0.0_wp)*dht*N*h*dht - (0.0_wp,2.0_wp)*e*P*ht
    end associate
  end subroutine

  pure subroutine halfmetal_interface_equation_a(this, a, g, gt, dg, dgt, r, rt)
    ! Calculate residuals from the boundary conditions at the left interface.
    class(halfmetal),          intent(in   ) :: this
    type(green),               intent(in   ) :: a
    type(spin),                intent(in   ) :: g, gt, dg, dgt
    type(spin),                intent(inout) :: r, rt

    ! Diagonal components: use spin-active boundary conditions
    if (associated(this%material_a)) then
      ! Interface is tunneling
      call this%interface_tunnel_a(a, g, gt, dg, dgt, r, rt)
      call this%interface_spinactive_a(a, g, gt, dg, dgt, r, rt)
    else
      ! Interface is vacuum
      call this%interface_vacuum_a(g, gt, dg, dgt, r, rt)
    end if

    ! Off-diagonal components: use boundary conditions g,gt=0
    r  % matrix(1,2)  = g  % matrix(1,2)
    r  % matrix(2,1)  = g  % matrix(2,1)
    rt % matrix(1,2)  = gt % matrix(1,2)
    rt % matrix(2,1)  = gt % matrix(2,1)
  end subroutine

  pure subroutine halfmetal_interface_equation_b(this, b, g, gt, dg, dgt, r, rt)
    ! Calculate residuals from the boundary conditions at the right interface.
    class(halfmetal),          intent(in   ) :: this
    type(green),               intent(in   ) :: b
    type(spin),                intent(in   ) :: g, gt, dg, dgt
    type(spin),                intent(inout) :: r, rt

    ! Diagonal components: use spin-active boundary conditions
    if (associated(this%material_b)) then
      ! Interface is tunneling
      call this%interface_tunnel_b(b, g, gt, dg, dgt, r, rt)
      call this%interface_spinactive_b(b, g, gt, dg, dgt, r, rt)
    else
      ! Interface is vacuum
      call this%interface_vacuum_b(g, gt, dg, dgt, r, rt)
    end if

    ! Off-diagonal components: use boundary conditions g,gt=0
    r  % matrix(1,2)  = g  % matrix(1,2)
    r  % matrix(2,1)  = g  % matrix(2,1)
    rt % matrix(1,2)  = gt % matrix(1,2)
    rt % matrix(2,1)  = gt % matrix(2,1)
  end subroutine

  impure subroutine halfmetal_update_prehook(this)
    ! Code to execute before running the update method of a class(halfmetal) object.
    class(halfmetal), intent(inout) :: this

    ! Verify that a polarization is defined
    if (abs(this % polarization) < eps) then
      call error('Tried to update a halfmetal with no defined polarization!')
    end if

    ! Update the polarization matrix
    this % P % matrix(1,1) = 2/(1 + eps + this%polarization)
    this % P % matrix(1,2) = 0
    this % P % matrix(2,1) = 0
    this % P % matrix(2,2) = 2/(1 + eps - this%polarization)

    ! Update the left  interface parameters
    this % magnetization_a = [0,0,1]
    this % polarization_a  = this % polarization

    ! Update the right interface parameters
    this % magnetization_b = [0,0,1]
    this % polarization_b  = this % polarization

    ! Call the superclass prehook
    call this%conductor%update_prehook

    ! Modify the type string
    this%type_string = color_red // 'HALFMETAL' // color_none
  end subroutine

  impure subroutine halfmetal_update_posthook(this)
    ! Code to execute after running the update method of a class(halfmetal) object.
    class(halfmetal), intent(inout) :: this
    real(wp)                        :: error
    integer                         :: n, m

    ! Call the superclass posthook
    call this%conductor%update_posthook

    ! Perform sanity checks
    error = 0
    do m=1,size(this%location)
      do n=1,size(this%energy)
        ! Find the maximum off-diagonal element
        error = max(abs(this%propagator(n,m)%g  % matrix(1,2)), error)
        error = max(abs(this%propagator(n,m)%g  % matrix(2,1)), error)
        error = max(abs(this%propagator(n,m)%gt % matrix(1,2)), error)
        error = max(abs(this%propagator(n,m)%gt % matrix(2,1)), error)
      end do
    end do

    ! Status information
    if (this%information >= 0 .and. .not. this % lock) then
      write(stdout,'(4x,a,f10.8,a)') 'Max error:  ',error,'                                        '
      flush(stdout)
    end if
  end subroutine

  !--------------------------------------------------------------------------------!
  !                      IMPLEMENTATION OF UTILITY METHODS                         !
  !--------------------------------------------------------------------------------!

  impure subroutine halfmetal_conf(this, key, val)
    !! Configure a material property based on a key-value pair.
    class(halfmetal), intent(inout) :: this
    character(*),     intent(in   ) :: key
    character(*),     intent(in   ) :: val
    real(wp)                        :: tmp

    select case(key)
      case ('polarization')
        read(val,*) this % polarization
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
      case default
        call material_conf(this, key, val)
    end select
  end subroutine
end module
