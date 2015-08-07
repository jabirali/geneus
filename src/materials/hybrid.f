! This module defines a set of subroutines and functions that are useful for working with multilayer hybrid structures.
! The procedures include the subroutine 'connect' for creating interfaces between class(material) objects; subroutines 
! 'init_all' and 'update_all' for manipulating the internal states of all materials in a hybrid structure; and a set of
! functions that are useful for initializing arguments that will be passed to the class(material) constructor methods.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-11
! Updated: 2015-07-29

module mod_hybrid
  use mod_system
  use mod_material
  use mod_conductor
  use mod_superconductor
  use mod_ferromagnet
  implicit none
contains

  !--------------------------------------------------------------------------------!
  !                PROCEDURES FOR ASSEMBLING MULTILAYER STRUCTURES                 !
  !--------------------------------------------------------------------------------!

  subroutine connect(material_a, material_b, conductance_a, conductance_b)
    ! This subroutine connects two class(material) materials by a tunneling interface, and may
    ! therefore be used to assemble individual material layers to a multilayer hybrid structure.
    class(material), target, intent(inout) :: material_a      ! This object represents the left  material
    class(material), target, intent(inout) :: material_b      ! This object represents the right material
    real(dp),                intent(in)    :: conductance_a   ! Tunneling conductance of the interface (relative to the left  bulk conductance)
    real(dp),                intent(in)    :: conductance_b   ! Tunneling conductance of the interface (relative to the right bulk conductance)
    
    ! Update the internal material pointers
    material_a % material_b => material_b
    material_b % material_a => material_a

    ! Update the interface parameters
    material_a % conductance_b = conductance_a
    material_b % conductance_a = conductance_b
  end subroutine



  !--------------------------------------------------------------------------------!
  !              PROCEDURES FOR CLASS(CONDUCTOR) CONSTRUCTOR ARGUMENTS             !
  !--------------------------------------------------------------------------------!

  pure subroutine energy_range(array, coupling, maximum, padding)
    ! Initializes an array of energies, which can be passed on to class(material) constructor methods.  The initialized
    ! values depend on the optional arguments.  If none of these arguments are provided, then the array is initilized to
    ! linearly spaced values in the range [0.0,1.5]. If the argument 'positive' is set to false, then the array includes
    ! negative values as well, resulting in linearly spaced values in the range [-1.5,1.5]. The default limit of 1.5 can
    ! be changed using the argument 'maximum'.  Finally, if the coupling constant 'coupling' is provided, the array will
    ! be padded with 100 linearly spaced energies up to the Debye cutoff cosh(1/coupling). If 'positive' is set to false,
    ! this becomes 100 positive and 100 negative energies. The parameter 'padding' can be used to change the default 100.
    !
    ! TODO:  the argument 'positive' has not been implemented yet,  and the other arguments are not checked for validity.

    real(dp), intent(out)          :: array(:)
    real(dp), optional, intent(in) :: maximum
    real(dp), optional, intent(in) :: coupling
    integer,  optional, intent(in) :: padding

    real(dp)                       :: maximum_
    integer                        :: padding_
    integer                        :: n

    ! Handle optional arguments
    if (present(maximum)) then
      maximum_ = maximum
    else
      maximum_ = 1.5_dp
    end if

    if (present(padding)) then
      padding_ = padding
    else
      padding_ = 100
    end if

    ! Initialize the energy array
    if (present(coupling)) then
      ! Positive energies from zero to 'maximum'
      do n = 1,size(array)-padding_
        array(n) = (n-1) * (maximum_/(size(array)-padding_))
      end do
      ! Positive energies from 'maximum' to the Debye cutoff
      do n = 1,padding_
        array(size(array)-padding_+n) = maximum_ + n * (cosh(1.0_dp/coupling)-maximum_)/padding_
      end do
    else
      ! Positive energies from zero to 'maximum'
      do n = 1,size(array)
        array(n) = (n-1) * (maximum_/(size(array)-1))
      end do
    end if
  end subroutine

  pure function exchange_xy(strength, angle) result(field)
    ! This function returns a vector that describes an exchange field in the xy-plane,
    ! where the input arguments describe the exchange field using polar coordinates.
    real(dp), intent(in) :: strength
    real(dp), intent(in) :: angle
    real(dp)             :: field(3)

    field(1) = strength * cos(angle)
    field(2) = strength * sin(angle)
    field(3) = 0.0_dp
  end function

  pure function spinorbit_xy(strength, angle, alpha, beta) result(field)
    ! This function returns an SU(2) vector that describes a Rashba--Dresselhaus coupling
    ! in the xy-plane. The coupling constants are given in polar coordinates, so that the
    ! Rashba constant is strength*sin(angle), and the Dresselhaus one strength*cos(angle).
    ! Alternatively, the Rashba constant alpha and Dresselhaus constant beta can also be
    ! explicitly specified.
    real(dp), intent(in), optional :: strength
    real(dp), intent(in), optional :: angle
    real(dp), intent(in), optional :: alpha
    real(dp), intent(in), optional :: beta
    type(spin)                     :: field(3)

    ! Initialize to zero
    field(:) = spin(0)

    ! Polar parametrization
    if (present(strength) .and. present(angle)) then
      field(1) = (+strength) * (cos(angle)*pauli1 + sin(angle)*pauli2)
      field(2) = (-strength) * (cos(angle)*pauli2 + sin(angle)*pauli1)
    end if

    ! Explicit Rashba coefficient
    if (present(alpha)) then
      field(1) = (+beta)*pauli1 + (+alpha)*pauli2
    end if

    ! Explicit Dresselhaus coefficient
    if (present(beta)) then
      field(2) = (-beta)*pauli2 + (-alpha)*pauli1
    end if
  end function
end module
