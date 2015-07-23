! CR: 2015-07-23
! UP   -- " --

module mod_multilayer
  use mod_system
  use mod_conductor
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
    do while (.not. associated(p % material_b, m))
      p => p % material_b
      call p % update
    end do

    ! Restart at the specified material
    p => m

    ! Move to the right until a vacuum interface is located, and update the states of all materials on the way
    do while (associated(p % material_b))
      p => p % material_b
      call p % update
    end do

    ! Move to the left until we get back to the start point, and update the states of all materials on the way
    do while (.not. associated(p % material_a, m))
      p => p % material_a
      call p % update
    end do

    ! Finally, update the state of the specified material itself
    call m % update
  end subroutine

  pure subroutine energy_range_positive(array, coupling)
    ! Initializes an energy array to values from zero to the Debye cutoff,
    ! where the Debye cutoff is calculated from the BCS coupling constant.
    real(dp), intent(out) :: array(:)  ! Array that will be initialized
    real(dp), intent(in)  :: coupling  ! BCS coupling constant
    integer               :: n         ! Loop variable

    ! Let the first N-100 elements span the energy range from zero to 1.5Δ
    do n=1,size(array)-100
      array(n) = (n-1)*(1.5_dp/(size(array)-100))
    end do

    ! Let the last 100 elements span the energy range from 1.5Δ to the cutoff
    do n=1,100
      array(size(array)-100+n) = 1.5_dp + n*(cosh(1.0_dp/coupling) - 1.5_dp)/100.0_dp
    end do
  end subroutine
end module
