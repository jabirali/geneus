!! This module defines a data type 'structure', which is useful for constructing and using multilayer hybrid structures.
!! @NOTE: This module only exports a single object 'structure', which encapsulates all other required objects/routines!
!! @TODO: This module is supposed to replace/supersede mod_multilayer, mod_hybrid, mod_option, and the older mod_config.
!!
!! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
!! Created: 2016-03-08
!! Updated: 2016-03-08

module mod_structure
  use mod_math
  use mod_stdio
  use mod_material
  use mod_conductor
  use mod_superconductor
  use mod_ferromagnet
  implicit none
  private

  public structure

  type :: structure
    class(material), pointer :: a => null()
    class(material), pointer :: b => null()
  contains
    procedure :: push_back => structure_push_back
    procedure :: conf_back => structure_conf_back
    procedure :: update    => structure_update
  end type
contains
  impure subroutine structure_push_back(struct, string)
    !! Constructs a new class(material) object at the bottom of the multilayer stack.
    class(structure), intent(inout) :: struct
    character(*),     intent(in   ) :: string

    if (.not. associated(struct % b)) then
      ! This is the first layer in the structure
      call material_allocate(struct % b, string)
      call material_construct(struct % b)
      struct % a => struct % b
    else
      ! This is not the first layer in the structure
      call material_allocate(struct % b % material_b, string)
      call material_construct(struct % b % material_b)
      struct % b % material_b % material_a => struct % b
      struct % b => struct % b % material_b
    end if
  contains
    impure subroutine material_allocate(ptr, str)
      !! Allocates memory for a new material layer.
      class(material), pointer, intent(inout) :: ptr
      character(*),             intent(in   ) :: str

      select case(str)
        case('ferromagnet')
          allocate(ferromagnet    :: ptr)
        case('superconductor')
          allocate(superconductor :: ptr)
        case('conductor')
          allocate(conductor      :: ptr)
        case default
          call error('Material type "' // str // '" unknown!')
      end select
    end subroutine

    impure subroutine material_construct(ptr)
      !! Constructs and initializes a new material layer.
      class(material), pointer, intent(inout) :: ptr

      select type(ptr)
        type is (ferromagnet)
          ptr = ferromagnet(30.0_wp, exchange = [10.0_wp,0.0_wp,0.0_wp])
        type is (superconductor)
          ptr = superconductor(30.0_wp)
        type is (conductor)
          ptr = conductor(30.0_wp)
        class default
          call error('Attempted to construct an unsupported material type!')
      end select
    end subroutine
  end subroutine

  impure subroutine structure_conf_back(struct, key, val)
    !! Configures the last material pushed to the multilayer stack.
    class(structure), intent(inout) :: struct
    character(*),     intent(in   ) :: key
    character(*),     intent(in   ) :: val

    if (associated(struct % b)) then
      call struct % b % conf(key, val)
    else
      call error('Attempted to configure a non-existant material!')
    end if
  end subroutine

  impure subroutine structure_update(struct)
    !! Updates the state of the entire multilayer stack.
    class(structure), target  :: struct
    class(material),  pointer :: ptr

    ! Initialize the material pointer to the top of the stack
    ptr => struct % a

    ! Update all materials (going down)
    do while (associated(ptr % material_b))
      ptr => ptr % material_b
      call ptr % update
    end do

    ! Update all materials (going up)
    do while (associated(ptr % material_a))
      ptr => ptr % material_a
      call ptr % update
    end do
  end subroutine
end module
