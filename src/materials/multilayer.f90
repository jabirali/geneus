! CR: 2015-07-23
! UP   -- " --

module mod_multilayer
  use mod_system
  use mod_conductor
  use mod_superconductor
  use mod_ferromagnet
  implicit none
contains
  subroutine connect(material_a, material_b, conductance_a, conductance_b)
    ! This subroutine connects two class(conductor) materials by a tunneling interface, and may
    ! therefore be used to assemble individual material layers to a multilayer hybrid structure.
    class(conductor), target, intent(inout) :: material_a      ! This object represents the left  material
    class(conductor), target, intent(inout) :: material_b      ! This object represents the right material
    real(dp),                 intent(in)    :: conductance_a   ! Tunneling conductance of the interface (relative to the left  bulk conductance)
    real(dp),                 intent(in)    :: conductance_b   ! Tunneling conductance of the interface (relative to the right bulk conductance)
    
    ! Update the internal material pointers
    material_a % material_b => material_b
    material_b % material_a => material_a

    ! Update the interface parameters
    material_a % conductance_b = conductance_a
    material_b % conductance_a = conductance_b
  end subroutine

  subroutine initialize_all(m, gap)
    ! This subroutine is used to initialize the physical state of a multilayer hybrid system to a BCS superconductor,
    ! where the system should have been assembled by previous calls to the routine 'connect'.  The subroutine takes a
    ! single class(conductor) object 'm' as its argument,  and initializes the entire associated multilayer structure
    ! to a BCS state with a given gap by jumping from layer to layer, and calling their respective initialize-methods.
    class(conductor), target     :: m
    class(conductor), pointer    :: p
    complex(dp),      intent(in) :: gap

    ! Initialize the specified material
    call m % initialize(gap)

    ! Start at the specified material, move left  until a vacuum interface is located, and initialize all layers on the way
    p => m
    do while (associated(p % material_a))
      p => p % material_a
      call p % initialize(gap)
    end do

    ! Start at the specified material, move right until a vacuum interface is located, and initialize all layers on the way
    p => m
    do while (associated(p % material_b))
      p => p % material_b
      call p % initialize(gap)
    end do
  end subroutine

  subroutine update_all(m)
    ! This subroutine is used to update the physical state of a multilayer hybrid system, where the system should have
    ! been assembled by previous calls to the routine 'connect'. The subroutine takes a single class(conductor) object
    ! 'm' as its argument,  and first updates the states of all layers to the left of 'm',  then updates the states of 
    ! all materials to the right of 'm', and then finally updates the state of 'm' itself. Since 'm' will be the final
    ! material updated, and therefore the material with the most accurately described state, this should be a layer of
    ! particular physical significance, such as the superconductor investigated in a critical temperature calculation.
    class(conductor), target  :: m
    class(conductor), pointer :: p

    ! Start at the specified material
    p => m

    ! Move to the left until a vacuum interface is located, and update the states of all materials on the way
    do while (associated(p % material_a))
      p => p % material_a
      call p % update
    end do

    ! Move to the right until we get back to the start point, and update the states of all materials on the way
    if (.not. associated(p, m)) then
      do while (associated(p % material_b) .and. .not. associated(p % material_b, m))
        p => p % material_b
        call p % update
      end do
    end if

    ! Restart at the specified material
    p => m

    ! Move to the right until a vacuum interface is located, and update the states of all materials on the way
    do while (associated(p % material_b))
      p => p % material_b
      call p % update
    end do

    ! Move to the left until we get back to the start point, and update the states of all materials on the way
    if (.not. associated(p, m)) then
      do while (associated(p % material_a) .and. .not. associated(p % material_a, m))
        p => p % material_a
        call p % update
      end do
    end if

    ! Finally, update the state of the specified material itself
    call m % update
  end subroutine

  pure subroutine energy_range(array, coupling, maximum, padding)
    ! Initializes an array of energies, which can be passed on to class(conductor) constructor methods.  The initialized
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
        array(n) = (n-1) * (maximum_/size(array))
      end do
    end if
  end subroutine
end module
