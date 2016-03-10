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
  use mod_halfmetal
  implicit none
  private

  type, public :: structure
    class(material), pointer :: a => null()
    class(material), pointer :: b => null()
  contains
    procedure :: push_back     => structure_push_back
    procedure :: conf_back     => structure_conf_back
    procedure :: map           => structure_map
    procedure :: init          => structure_init
    procedure :: save          => structure_save
    procedure :: load          => structure_load
    procedure :: update        => structure_update
    procedure :: difference    => structure_difference
    procedure :: temperature   => structure_temperature
    procedure :: write_density => structure_write_density
    procedure :: write_current => structure_write_current
    procedure :: write_gap     => structure_write_gap
  end type

  interface structure
    module procedure structure_construct
  end interface
contains
  impure subroutine structure_push_back(this, string)
    !! Constructs a new class(material) object at the bottom of the multilayer stack.
    class(structure), intent(inout) :: this
    character(*),     intent(in   ) :: string

    if (.not. associated(this % b)) then
      ! This is the first layer in the structure
      call material_allocate(this % b, string)
      call material_construct(this % b)
      this % a => this % b
    else
      ! This is not the first layer in the structure
      call material_allocate(this % b % material_b, string)
      call material_construct(this % b % material_b)
      this % b % material_b % material_a => this % b
      this % b => this % b % material_b
    end if
  contains
    impure subroutine material_allocate(ptr, str)
      !! Allocates memory for a new material layer.
      class(material), pointer, intent(inout) :: ptr
      character(*),             intent(in   ) :: str

      select case(str)
        case('halfmetal')
          allocate(halfmetal      :: ptr)
        case('ferromagnet')
          allocate(ferromagnet    :: ptr)
        case('superconductor')
          allocate(superconductor :: ptr)
        case('conductor')
          allocate(conductor      :: ptr)
        case default
          call error('Material type "' // trim(str) // '" unknown!')
      end select
    end subroutine

    impure subroutine material_construct(ptr)
      !! Constructs and initializes a new material layer.
      class(material), pointer, intent(inout) :: ptr

      select type(ptr)
        type is (halfmetal)
          ptr = halfmetal()
        type is (ferromagnet)
          ptr = ferromagnet()
        type is (superconductor)
          ptr = superconductor()
        type is (conductor)
          ptr = conductor()
        class default
          call error('Attempted to construct an unsupported material type!')
      end select
    end subroutine
  end subroutine

  impure subroutine structure_conf_back(this, key, val)
    !! Configures the last material pushed to the multilayer stack.
    class(structure), intent(inout) :: this
    character(*),     intent(in   ) :: key
    character(*),     intent(in   ) :: val

    if (associated(this % b)) then
      call this % b % conf(key, val)
    else
      call error('Attempted to configure a non-existant material!')
    end if
  end subroutine

  impure subroutine structure_map(this, routine, loop)
    !! Maps a subroutine onto each element of the multilayer stack.
    class(structure), target  :: this
    class(material),  pointer :: ptr
    external                  :: routine
    logical, optional         :: loop

    ! Initialize the material pointer to the top of the stack
    ptr => this % a

    if (.not. associated(ptr % material_b)) then
      ! Operate on this layer if it is the only one
      call routine(ptr)
    else if (.not. present(loop)) then
      ! Iterate through all layers once if 'loop' is not set
      do while (associated(ptr % material_b))
        ptr => ptr % material_b
        call routine(ptr)
      end do
    else
      ! Iterate through all layers but the first (going down)
      do while (associated(ptr % material_b))
        ptr => ptr % material_b
        call routine(ptr)
      end do
      ! Iterate through all layers but the last  (going up)
      do while (associated(ptr % material_a))
        ptr => ptr % material_a
        call routine(ptr)
      end do
    end if
  end subroutine

  impure subroutine structure_init(this, gap)
    !! Initializes the state of the entire multilayer stack.
    class(structure), target  :: this
    complex(wp)               :: gap

    ! Initialize all material states
    call this % map(init)
  contains
    subroutine init(m)
      class(material) :: m
      call m % init(gap)
    end subroutine
  end subroutine

  impure subroutine structure_save(this)
    !! Saves the state of the entire multilayer stack.
    class(structure), target  :: this

    ! Save all material states
    call this % map(save)
  contains
    subroutine save(m)
      class(material) :: m
      call m % save
    end subroutine
  end subroutine

  impure subroutine structure_load(this)
    !! Loads the saved state of the multilayer stack.
    class(structure), target  :: this

    ! Load all material states
    call this % map(load)
  contains
    subroutine load(m)
      class(material) :: m
      call m % load
    end subroutine
  end subroutine

  impure subroutine structure_update(this, freeze)
    !! Updates the state of the entire multilayer stack.
    class(structure), target   :: this
    logical,          optional :: freeze

    ! Update all materials in a loop-pattern
    call this % map(update, loop = .true.)
  contains
    subroutine update(m)
      class(material) :: m
      call m % update(freeze)
    end subroutine
  end subroutine

  impure function structure_difference(this) result(difference)
    !! Checks how much the multilayer stack has changed recently.
    class(structure), target  :: this
    real(wp)                  :: difference

    ! Check all materials
    difference = 0
    call this % map(check)
  contains
    subroutine check(m)
      class(material) :: m
      difference = max(difference, m % difference)
    end subroutine
  end function

  impure subroutine structure_temperature(this, temperature)
    !! Modifies the temperature of the multilayer stack.
    class(structure), target  :: this
    real(wp)                  :: temperature

    ! Update the temperature of all materials
    call this % map(set)
  contains
    subroutine set(m)
      class(material) :: m
      m % temperature = temperature
    end subroutine
  end subroutine


  impure subroutine structure_write_density(this, file)
    ! Writes the density of states as a function of position and energy to a given output file.
    class(structure), target, intent(in) :: this
    character(*),             intent(in) :: file
    integer                              :: unit
    integer                              :: iostat
    class(material),  pointer            :: ptr
    real(wp)                             :: a, b, x
    integer                              :: n, m

    ! Open output file
    open(newunit = unit, file = file, iostat = iostat, action = 'write', status = 'replace')
    if (iostat /= 0) then
      call error('failed to open output file "' // file // '"!')
    end if

    ! Initialize variables
    ptr => this % a
    b = 0

    ! Traverse all materials (going down)
    do while (associated(ptr))
      ! Calculate the endpoints of this layer
      a = b
      b = b + 1/sqrt(ptr % thouless)

      ! Write out the density of states in this layer
      do m=1,size(ptr % location)
        x = a+sqrt(eps) + ((b-sqrt(eps))-(a+sqrt(eps))) * ptr % location(m)
        do n=size(ptr % energy),1,-1
          ! Negative energies
          write(unit,*) x, -ptr % energy(n), ptr % density(n,m)
        end do
        do n=1,size(ptr % energy),+1
          ! Positive energies
          write(unit,*) x, +ptr % energy(n), ptr % density(n,m)
        end do
      end do

      ! Move on to the next material layer
      ptr => ptr % material_b
    end do

    ! Close output file
    close(unit = unit)
  end subroutine

  impure subroutine structure_write_current(this, file)
    ! Writes the charge and spin currents as a function of position to a given output file.
    class(structure), target, intent(in) :: this
    character(*),             intent(in) :: file
    integer                              :: unit
    integer                              :: iostat
    class(material),  pointer            :: ptr
    real(wp)                             :: a, b, x
    integer                              :: m

    ! Open output file
    open(newunit = unit, file = file, iostat = iostat, action = 'write', status = 'replace')
    if (iostat /= 0) then
      call error('failed to open output file "' // file // '"!')
    end if

    ! Initialize variables
    ptr => this % a
    b = 0

    ! Traverse all materials (going down)
    do while (associated(ptr))
      ! Calculate the endpoints of this layer
      a = b
      b = b + 1/sqrt(ptr % thouless)

      ! Write out the currents in this layer
      do m=1,size(ptr % location)
        x = a+sqrt(eps) + ((b-sqrt(eps))-(a+sqrt(eps))) * ptr % location(m)
        write(unit,*) x, ptr % current(:,m)
      end do

      ! Move on to the next material layer
      ptr => ptr % material_b
    end do

    ! Close output file
    close(unit = unit)
  end subroutine

  impure subroutine structure_write_gap(this, file)
    ! Writes the superconducting order parameter as a function of position to a given output file.
    class(structure), target, intent(in) :: this
    character(*),             intent(in) :: file
    integer                              :: unit
    integer                              :: iostat
    class(material),  pointer            :: ptr
    real(wp)                             :: a, b, x
    integer                              :: m

    ! Open output file
    open(newunit = unit, file = file, iostat = iostat, action = 'write', status = 'replace')
    if (iostat /= 0) then
      call error('failed to open output file "' // file // '"!')
    end if

    ! Initialize variables
    ptr => this % a
    b = 0

    ! Traverse all materials (going down)
    do while (associated(ptr))
      ! Calculate the endpoints of this layer
      a = b
      b = b + 1/sqrt(ptr % thouless)

      ! Write out the gap in this layer
      do m=1,size(ptr % location)
        x = a+sqrt(eps) + ((b-sqrt(eps))-(a+sqrt(eps))) * ptr % location(m)
        select type (ptr)
          class is (superconductor)
            write(unit,*) x, abs(ptr%gap(m)), atan2(im(ptr%gap(m)),re(ptr%gap(m)))/pi
          class default
            write(unit,*) x, 0.0_wp, 0.0_wp
        end select
      end do

      ! Move on to the next material layer
      ptr => ptr % material_b
    end do

    ! Close output file
    close(unit = unit)
  end subroutine

  impure function structure_construct(file) result(this)
    !! Constructs a multilayer stack from a configuration file.
    type(structure)          :: this
    character(*), intent(in) :: file
    integer                  :: unit
    integer                  :: iostat
    character(len=132)       :: str, lhs, rhs, arg
    integer                  :: line, i, j, k

    ! Initialize variables
    line   = 0
    unit   = 0
    iostat = 0

    ! Open the config file
    open(newunit = unit, file = file, iostat = iostat, action = 'read', status = 'old')
    if (iostat /= 0) then
      call error('failed to open configuration file "' // file // '"!')
    end if

    ! Read the config file
    do while (iostat == 0)
      ! Increase line counter
      line = line + 1

      ! Read one line from file
      str = ''
      read(unit, '(a)', iostat = iostat) str

      ! Strip comments from line
      i = scan(str, '#')
      if ( i > 0 ) then
        str = str(:i-1)
      end if

      ! Strip whitespace from line
      str = trim(str)

      ! Ignore simulation properties
      ! @TODO: implementation
      if (scan(str, '@') /= 0) then
        cycle
      end if
 
      ! Substitute command line arguments
      i = scan(str, '{')
      j = scan(str, '}')
      do while ( i > 0 .and. j-1 >= i+1 )
        read(str(i+1:j-1),*) k
        call get_command_argument(k, arg, status=k)
        if ( k /= 0 ) then
          call error('missing command line arguments.')
        end if
        str = str(:i-1) // trim(arg) // str(j+1:)
        i = scan(str, '{')
        j = scan(str, '}')
      end do

      ! Configure or construct a material
      i = scan(str, ':')
      if ( i > 0 ) then
        lhs = adjustl(str(:i-1))
        rhs = adjustl(str(i+1:))
        if (lhs /= '') then
          if (rhs /= '') then
            call this % conf_back(lhs, rhs)
          else
            call this % push_back(lhs)
          end if
          cycle
        end if
      end if

      ! Check if this line contains garbage
      if ( str /= '' ) then
        write(str,'(i0)') line
        call error('failed to parse line ' // trim(str) // ' in the config file.')
      end if
    end do

    ! Close the config file
    close(unit = unit)
  end function
end module
