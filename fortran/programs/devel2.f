module mod_structure
  use mod_hybrid
  use mod_stdio

  type :: structure
    class(material), pointer :: a => null()
    class(material), pointer :: b => null()
  contains
    procedure :: push   => structure_push
    procedure :: update => structure_update
  end type
contains
  impure subroutine structure_push(struct, string)
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

  impure subroutine structure_update(s)
    ! This subroutine updates the state of the entire hybrid structure.
    class(structure), target  :: s
    class(material),  pointer :: p

    ! Initialize the material pointer to the top of the stack
    p => s % a

    ! Update all materials (going down)
    do while (associated(p % material_b))
      p => p % material_b
      call p % update
    end do

    ! Update all materials (going up)
    do while (associated(p % material_a))
      p => p % material_a
      call p % update
    end do
  end subroutine
end module

program test
  use mod_hybrid
  use mod_material
  use mod_structure

  type(structure) :: bilayer

  call bilayer % push('superconductor')
  call bilayer % push('ferromagnet')

  call bilayer % a % conf('temperature', '0.10')
  call bilayer % a % conf('scattering',  '0.05')
  call bilayer % a % conf('length',      '0.75')
  call bilayer % a % conf('coupling',    '0.25')

  call bilayer % b % conf('temperature', '0.10')
  call bilayer % b % conf('scattering',  '0.05')
  call bilayer % b % conf('length',      '0.50')

  call bilayer % update
end program


!module mod_neoconf
!  use mod_material
!
!  type :: dict
!    !! Stores simple key-value mappings in a single-linked list.
!    character(len=132)   :: key  =  ''      ! Dictionary key
!    character(len=132)   :: val  =  ''      ! Dictionary value
!    class(dict), pointer :: next => null()  ! Next element
!  end type
!
!  type :: conf
!    !! Stores a multilayer configuration as a double-linked list.
!    character(len=132)   :: node =  ''      ! Name of the node
!    class(conf), pointer :: next => null()  ! Next node
!    class(conf), pointer :: prev => null()  ! Previous node
!    class(dict), pointer :: dict => null()  ! Dictionary
!  end type
!contains
!  impure subroutine neoconf_read()
!    !! Read one record from a config file, filtering whitespace and comments in the process, and respecting line continuations.
!    continue
!  end subroutine
!
!  impure subroutine neoconf_parse()
!    !! Parse one record from a config file, creating either a node or a dict mapping in the process.
!    continue
!  end subroutine
!
!  impure subroutine neoconf()
!    !! Perform the entire process.
!    continue
!  end subroutine
!end module

