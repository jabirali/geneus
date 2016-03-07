module mod_structure
  use mod_hybrid

  type :: structure
    class(material), allocatable :: top
  contains
    procedure :: update => structure_update
  end type
contains
  impure subroutine structure_update(s)
    ! This subroutine updates the state of the entire hybrid structure.
    class(structure), target  :: s
    class(material),  pointer :: p

    ! Initialize the material pointer to the top of the stack
    p => s % top

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
  use mod_neoconf

  type(superconductor) :: s
  type(ferromagnet)    :: f

  s = superconductor(30.0_wp)
  f = ferromagnet(30.0_wp)
  call connect(s,f)

  call f%conf('temperature', '0.10')
  call f%conf('scattering',  '0.05')
  call f%conf('length',      '0.50')

  call s%conf('temperature', '0.10')
  call s%conf('scattering',  '0.05')
  call s%conf('length',      '0.75')
  call s%conf('coupling',    '0.25')

  write(*,*) f%thouless
  write(*,*) f%temperature
  write(*,*) f%scattering

  call f%update
  call s%update
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

