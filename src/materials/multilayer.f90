
module mod_multilayer
  use mod_system
  use mod_conductor
  implicit none

  ! Connecting two class(conductor) objects by creating a tunneling interface between them
  interface connect
    module procedure connect_tunneling
  end interface
contains
  pure subroutine energy_range_positive(array, coupling)
    ! Initializes an energy array to values from zero to the Debye cutoff, where the cutoff is calculated from a BCS coupling constant.
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

  subroutine connect_tunneling(material_a, material_b, conductance_a, conductance_b)
    ! Creates a tunneling interface between two conductive materials.
    class(conductor), target, intent(inout) :: material_a      ! This object represents the left  material
    class(conductor), target, intent(inout) :: material_b      ! This object represents the right material
    real(dp),                 intent(in)    :: conductance_a   ! Tunneling conductance of the interface (relative to the left  bulk conductance)
    real(dp),                 intent(in)    :: conductance_b   ! Tunneling conductance of the interface (relative to the right bulk conductance)
    
    ! Update the internal material pointers
    material_a%material_b => material_b
    material_b%material_a => material_a

    ! Update the interface parameters
    material_a%conductance_b = conductance_a
    material_b%conductance_a = conductance_b
  end subroutine


  subroutine update(m)
    ! This subroutine is used to update the physical state of a multilayer hybrid system, where the hybrid system must
    ! have previously been constructed using the 'connect' subroutine. 
    class(conductor), target  :: m
    class(conductor), pointer :: p

    ! Start at the given material
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

    ! Restart at the given material
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

    ! Finally, update the state of the provided material
    call m % update
  end subroutine
end module
