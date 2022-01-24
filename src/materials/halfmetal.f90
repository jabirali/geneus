!> Author:   Jabir Ali Ouassou
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
  use :: stdio_m
  use :: condmat_m
  use :: material_m
  use :: conductor_m
  private

  ! Type declarations
  type, public, extends(conductor) :: halfmetal
    real(wp)            :: polarization = 0.0_wp                                  !! Spin-polarization of the ferromagnet
    type(spin), private :: P                                                      !! Polarization matrix
  contains
    procedure           :: conf                 => halfmetal_conf                 !! Configures the material parameters
    procedure           :: diffusion_equation   => halfmetal_diffusion_equation   !! Defines the Usadel diffusion equation
    procedure           :: diffusion_equation_a => halfmetal_diffusion_equation_a !! Boundary condition at the left  interface
    procedure           :: diffusion_equation_b => halfmetal_diffusion_equation_b !! Boundary condition at the right interface
    procedure           :: update_prehook       => halfmetal_update_prehook       !! Code to execute before calculating the propagators
    procedure           :: update_posthook      => halfmetal_update_posthook      !! Code to execute after  calculating the propagators
    procedure           :: update_density       => halfmetal_update_density       !! Calculates the density of states
  end type
contains

  !--------------------------------------------------------------------------------!
  !                     IMPLEMENTATION OF HALFMETAL EQUATIONS                      !
  !--------------------------------------------------------------------------------!

  pure subroutine halfmetal_diffusion_equation(this, p, e, z)
    !! Use the diffusion equation to calculate the second-derivatives
    !! of the Riccati parameters at an energy e and a position z.
    class(halfmetal), intent(in)    :: this
    complex(wp),      intent(in)    :: e
    real(wp),         intent(in)    :: z
    type(propagator), intent(inout) :: p
    type(spin)                      :: h, ht, dh, dht
    type(spin)                      :: N, Nt

    associate(  g => p % g,     gt => p % gt,  &
               dg => p % dg,   dgt => p % dgt, &
              d2g => p % d2g, d2gt => p % d2gt )

      ! Ensure that the Riccati parameters are diagonal
      h   = g   % matrix * pauli0 % matrix
      ht  = gt  % matrix * pauli0 % matrix
      dh  = dg  % matrix * pauli0 % matrix
      dht = dgt % matrix * pauli0 % matrix

      ! Calculate the normalization matrices
      N   = inverse( pauli0 - h*ht )
      Nt  = inverse( pauli0 - ht*h )

      ! Calculate the second-derivatives of the Riccati parameters
      associate(P => this % P)
        d2g  = (-2.0_wp,0.0_wp)*dh*Nt*ht*dh - (0.0_wp,2.0_wp)*e*P*h
        d2gt = (-2.0_wp,0.0_wp)*dht*N*h*dht - (0.0_wp,2.0_wp)*e*P*ht
      end associate
    end associate
  end subroutine

  pure subroutine halfmetal_diffusion_equation_a(this, p, a, r, rt)
    !! Calculate residuals from the boundary conditions at the left interface.
    class(halfmetal), intent(in)    :: this
    type(propagator), intent(in)    :: p, a
    type(spin),       intent(inout) :: r, rt

    ! Diagonal components: use regular spin-active boundary conditions
    call this % conductor % diffusion_equation_a(p, a, r, rt)

    ! Off-diagonal components: use boundary conditions g,gt=0
    r  % matrix(1,2)  = p % g  % matrix(1,2)
    r  % matrix(2,1)  = p % g  % matrix(2,1)
    rt % matrix(1,2)  = p % gt % matrix(1,2)
    rt % matrix(2,1)  = p % gt % matrix(2,1)
  end subroutine

  pure subroutine halfmetal_diffusion_equation_b(this, p, b, r, rt)
    !! Calculate residuals from the boundary conditions at the right interface.
    class(halfmetal), intent(in)    :: this
    type(propagator), intent(in)    :: p, b
    type(spin),       intent(inout) :: r, rt

    ! Diagonal components: use regular spin-active boundary conditions
    call this % conductor % diffusion_equation_b(p, b, r, rt)

    ! Off-diagonal components: use boundary conditions g,gt=0
    r  % matrix(1,2)  = p % g  % matrix(1,2)
    r  % matrix(2,1)  = p % g  % matrix(2,1)
    rt % matrix(1,2)  = p % gt % matrix(1,2)
    rt % matrix(2,1)  = p % gt % matrix(2,1)
  end subroutine

  subroutine halfmetal_update_prehook(this)
    !! Code to execute before running the update method of a class(halfmetal) object.
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
    this % spinactive_a % magnetization = [0,0,1]
    this % spinactive_a % polarization  = this % polarization

    ! Update the right interface parameters
    this % spinactive_b % magnetization = [0,0,1]
    this % spinactive_b % polarization  = this % polarization

    ! Call the superclass prehook
    call this%conductor%update_prehook

    ! Modify the type string
    this%type_string = color_red // 'HALFMETAL' // color_none
  end subroutine

  subroutine halfmetal_update_posthook(this)
    !! Code to execute after running the update method of a class(halfmetal) object.
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

    ! Update density of states
    call this%update_density

    ! Status information
    if (this%information >= 0 .and. this % order > 0) then
      write(stdout,'(6x,a,f10.8,a)') 'Max error:  ',error,'                                        '
    end if
  end subroutine

  pure subroutine halfmetal_update_density(this)
    !! Calculate the density of states in the halfmetal.
    class(halfmetal), intent(inout) :: this
    integer                         :: n, m

    ! Allocate memory if necessary
    if (.not. allocated(this%density)) then
      allocate(this%density(size(this%energy),size(this%location),0:7))
    end if

    ! Placeholder code
    this % density = inf
   
    ! ! Calculate the density of states at each position and energy
    ! ! TODO: Generalize this to work for strong ferromagnets too.
    ! if (this % polarization > 0) then
    !   do m=1,size(this%location)
    !     do n=1,size(this%energy)
    !       this % density(n,m) = 2 * re(this % propagator(n,m) % N % matrix(1,1)) - 1
    !     end do
    !   end do
    ! else 
    !   do m=1,size(this%location)
    !     do n=1,size(this%energy)
    !       this % density(n,m) = 2 * re(this % propagator(n,m) % N % matrix(2,2)) - 1
    !     end do
    !   end do
    ! end if
  end subroutine


  !--------------------------------------------------------------------------------!
  !                      IMPLEMENTATION OF UTILITY METHODS                         !
  !--------------------------------------------------------------------------------!

  subroutine halfmetal_conf(this, key, val)
    !! Configure a material property based on a key-value pair.
    use :: evaluate_m

    class(halfmetal), intent(inout) :: this
    character(*),     intent(in)    :: key
    character(*),     intent(in)    :: val
    real(wp)                        :: tmp

    select case(key)
      case ('polarization')
        call evaluate(val, this % polarization)
      case default
        call this % conductor % conf(key, val)
    end select
  end subroutine
end module
