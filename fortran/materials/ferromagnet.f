!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This module defines the data type 'ferromagnet', which models the physical state of a ferromagnet. The type is a
!> member of class(conductor), and thus inherits the internal structure and generic methods defined in conductor_m.

module ferromagnet_m
  use :: stdio_m
  use :: condmat_m
  use :: conductor_m
  private

  ! Type declaration
  type, public, extends(conductor)   :: ferromagnet
    real(wp),   allocatable          :: magnetization(:,:)                                    !! Magnetic exchange field as a function of position
    type(spin), allocatable, private :: h(:)                                                  !! Used by internal subroutines to handle exchange fields
  contains
    ! These methods define the class(material) interface
    procedure                        :: update_prehook     => ferromagnet_update_prehook      !! Code to execute before calculating the propagators

    ! These methods contain the equations that describe ferromagnets
    procedure                        :: diffusion_equation => ferromagnet_diffusion_equation  !! Diffusion equation
    procedure                        :: kinetic_equation   => ferromagnet_kinetic_equation    !! Kinetic equation

    ! These methods define miscellaneous utility functions
    procedure                        :: conf               => ferromagnet_conf                !! Configures material parameters
  end type
contains

  !--------------------------------------------------------------------------------!
  !                    IMPLEMENTATION OF FERROMAGNET METHODS                       !
  !--------------------------------------------------------------------------------!

  pure subroutine ferromagnet_diffusion_equation(this, p, e, z)
    !! Use the diffusion equation to calculate the second derivatives of the Riccati parameters at point z.
    class(ferromagnet), intent(in)    :: this
    complex(wp),        intent(in)    :: e
    real(wp),           intent(in)    :: z
    type(propagator),   intent(inout) :: p
    type(spin)                        :: h, ht
    real(wp)                          :: d
    integer                           :: n, m

    associate(  g => p % g,     gt => p % gt,  &
               dg => p % dg,   dgt => p % dgt, &
              d2g => p % d2g, d2gt => p % d2gt )


      ! Calculate the second derivatives of the Riccati parameters (conductor terms)
      call this % conductor % diffusion_equation(p, e, z)

      if (allocated(this%h)) then
        ! Calculate the index corresponding to the given location
        m = size(this%h) - 1          ! Number of array intervals
        n = floor(z*m + 1)            ! Nearest position in array

        ! Extract the exchange field terms at that point
        if (n <= 1) then
          ! Left edge of the material
          h  = this%h(1)
        else
          ! Linear interpolation from known values. The relative displacement d is defined
          ! as [z - location(n-1)]/[location(n) - location(n-1)], but assuming location(:)
          ! is a uniform array of values in the range [0,1], the below will be equivalent.
          d  = z*m - (n-2)
          h  = this%h(n-1)  + (this%h(n)  - this%h(n-1))  * d
        end if

        ! Find the corresponding tilde-conjugate
        ht = conjg(h)

        ! Calculate the second derivatives of the Riccati parameters (ferromagnet terms)
        associate( i => (0.00_wp,1.00_wp) )
          d2g  = d2g  - i * ( h  * g  - g  * ht )
          d2gt = d2gt + i * ( ht * gt - gt * h  )
        end associate
      end if

    end associate
  end subroutine

  pure subroutine ferromagnet_kinetic_equation(this, Gp, R, z)
    !! Calculate the self-energies in the kinetic equation.
    class(ferromagnet),           intent(in)    :: this
    type(propagator),             intent(in)    :: Gp
    real(wp), dimension(0:7,0:7), intent(inout) :: R
    real(wp),                     intent(in)    :: z
    type(nambu)                                 :: S
    real(wp)                                    :: d
    type(spin)                                  :: h, ht
    integer                                     :: n, m

    ! Call the superclass kinetic equation
    call this % conductor % kinetic_equation(Gp, R, z)

    if (allocated(this%h)) then
      ! Calculate the index corresponding to the given location
      m = size(this%h) - 1          ! Number of array intervals
      n = floor(z*m + 1)            ! Nearest position in array

      ! Extract the exchange field terms at that point
      if (n <= 1) then
        ! Left edge of the material
        h  = this%h(1)
      else
        ! Linear interpolation from known values. The relative displacement d is defined
        ! as [z - location(n-1)]/[location(n) - location(n-1)], but assuming location(:)
        ! is a uniform array of values in the range [0,1], the below will be equivalent.
        d  = z*m - (n-2)
        h  = this%h(n-1)  + (this%h(n)  - this%h(n-1))  * d
      end if

      ! Find the corresponding tilde-conjugate
      ht = conjg(h)

      ! Construct the self-energy matrix
      S % matrix(1:2,1:2) = h
      S % matrix(3:4,3:4) = ht

      ! Calculate the self-energy contribution
      R = R + re(Gp % selfenergy1(S))
    end if
  end subroutine

  impure subroutine ferromagnet_update_prehook(this)
    !! Updates the exchange field terms in the diffusion equation.
    class(ferromagnet), intent(inout) :: this ! Ferromagnet object that will be updated
    integer                           :: n    ! Loop variable

    ! Call the superclass prehook
    call this%conductor%update_prehook

    ! Rename the internal variables
    if (allocated(this%magnetization)) then
      ! Allocate space for internal variables
      if (.not. allocated(this%h)) then
        allocate(this%h (size(this%magnetization,2)))
      end if

      ! Update the internal variables
      associate(M => this % magnetization, h => this % h)
        do n = 1,size(M,2)
          h(n)  = (M(1,n)*pauli1 + M(2,n)*pauli2 + M(3,n)*pauli3)/(this%thouless)
        end do
      end associate
    end if

    ! Modify the type string
    if (allocated(this%magnetization)) then
      this%type_string = color_red // 'FERROMAGNET' // color_none
    end if
  end subroutine

  !--------------------------------------------------------------------------------!
  !                      IMPLEMENTATION OF UTILITY METHODS                         !
  !--------------------------------------------------------------------------------!

  impure subroutine ferromagnet_conf(this, key, val)
    !! Configure a material property based on a key-value pair.
    use :: evaluate_m

    class(ferromagnet), intent(inout)   :: this
    character(*),       intent(in)      :: key
    character(*),       intent(in)      :: val

    select case(key)
      case ('magnetization')
        block
          ! Allocate memory for the location array
          real(wp), allocatable, dimension(:) :: location
          allocate(location(1024 * size(this%location)))

          ! Discretize the locations in the material
          call linspace(location, 0.0_wp, 1.0_wp)

          ! Initialize the magnetic exchange field
          call evaluate(val, location, this % magnetization)

          ! Deallocate the field if it is negligible
          if (norm2(this % magnetization) < sqrt(eps)) then
            deallocate(this % magnetization)
          end if

          ! Deallocate the location array
          deallocate(location)
        end block
      case default
        ! Pass this option to the superclass
        call this % conductor % conf(key, val)
    end select
  end subroutine
end module
