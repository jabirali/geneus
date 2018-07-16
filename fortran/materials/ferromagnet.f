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
    real(wp),   allocatable          :: zeeman                                                !! How easy the material is magnetized by spin accumulation
    real(wp),   allocatable          :: exchange(:,:)                                         !! Magnetic exchange field as a function of position
    real(wp),   allocatable, private :: z(:)                                                  !! Used by internal subroutines to handle location data
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
    class(ferromagnet),              intent(in)    :: this
    type(propagator),                intent(in)    :: Gp
    complex(wp), dimension(0:7,0:7), intent(inout) :: R
    real(wp),                        intent(in)    :: z
    type(nambu)                                    :: S
    real(wp)                                       :: d
    type(spin)                                     :: h, ht
    integer                                        :: n, m

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
      R = R + Gp % selfenergy1(S)
    end if
  end subroutine

  impure subroutine ferromagnet_update_prehook(this)
    !! Updates the exchange field terms in the diffusion equation.
    class(ferromagnet), intent(inout)     :: this          ! Ferromagnet object that will be updated
    real(wp), allocatable, dimension(:,:) :: magnetization ! Locally calculated effective magnetization
    real(wp), allocatable, dimension(:,:) :: interpolation ! High-resolution magnetization interpolation
    integer                               :: n             ! Loop variable

    ! Call the superclass prehook
    call this % conductor % update_prehook

    ! Update magnetization matrices
    if (allocated(this % exchange) .or. allocated(this % zeeman)) then
      ! Allocate and initialize workspace
      allocate(magnetization(3, size(this % location)))
      allocate(interpolation(3, size(this % location) * 4096))
      if (.not. allocated(this % h)) then
        allocate(this % h(size(interpolation, 2)))
        allocate(this % z(size(interpolation, 2)))
        call linspace(this % z, this % location(1), this % location(size(this % location)))
      end if

      ! Calculate the effective magnetization
      magnetization = 0
      if (allocated(this % exchange)) then
        magnetization = magnetization + (this % exchange(1:3,:))
      end if
      if (allocated(this % zeeman)) then
        magnetization = magnetization + (this % accumulation(1:3,:)) * (this % zeeman)
      end if

      ! High-resolution interpolation
      interpolation(1,:) = interpolate(this % location, magnetization(1,:), this % z)
      interpolation(2,:) = interpolate(this % location, magnetization(2,:), this % z)
      interpolation(3,:) = interpolate(this % location, magnetization(3,:), this % z)

      ! Update the internal variables
      do n = 1,size(interpolation,2)
        write(*,*) interpolation(1:3,n)
        this % h(n) = (interpolation(1,n)*pauli1 + interpolation(2,n)*pauli2 + interpolation(3,n)*pauli3)/(this % thouless)
      end do

      ! Clean up
      deallocate(magnetization, interpolation)
    end if

    ! Modify the type string
    if (allocated(this % exchange) .or. allocated(this % zeeman)) then
      this%type_string = color_red // 'MAGNET' // color_none
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
        call evaluate(val, this % location, this % exchange)
      case ('zeeman')
        allocate(this % zeeman)
        call evaluate(val, this % zeeman)
      case default
        ! Pass this option to the superclass
        call this % conductor % conf(key, val)
    end select
  end subroutine
end module
