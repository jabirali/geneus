!! This module defines a data type 'structure', which is useful for constructing and using multilayer hybrid structures.
!! @NOTE: This module only exports a single object 'structure', which encapsulates all other required objects/routines!
!!
!! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
!! Created: 2016-03-08
!! Updated: 2016-03-23

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

  public conductor, superconductor, ferromagnet, halfmetal

  type, public :: structure
    class(material), pointer :: a => null()
    class(material), pointer :: b => null()
  contains
    procedure :: push            => structure_push
    procedure :: conf            => structure_conf
    procedure :: map             => structure_map
    procedure :: init            => structure_init
    procedure :: save            => structure_save
    procedure :: load            => structure_load
    procedure :: update          => structure_update
    procedure :: count           => structure_count
    procedure :: difference      => structure_difference
    procedure :: temperature     => structure_temperature
    procedure :: write_density   => structure_write_density
    procedure :: write_current   => structure_write_current
    procedure :: write_gap       => structure_write_gap
  end type

  interface structure
    module procedure structure_construct
  end interface
contains
  impure subroutine structure_push(this, string)
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

  impure subroutine structure_conf(this, key, val)
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
    logical                   :: loop_
    integer                   :: n

    ! Count the number of materials
    n = 0

    ! Determine whether or not we need to loop
    loop_ = .false.
    if (present(loop)) then
      loop_ = loop
    end if

    ! Determine where to start the traversation
    call top(ptr)
    if (loop_) then
      call next(ptr)
      n = n + 1
    end if

    ! Traverse the structure from top to bottom
    do while (associated(ptr))
      call routine(ptr)
      call next(ptr)
      n = n + 1
    end do

    ! Traverse the structure from bottom to top
    if (loop_) then
      call bottom(ptr)
      if (n > 1) then
        call prev(ptr)
      end if
      do while (associated(ptr))
        call routine(ptr)
        call prev(ptr)
      end do
    end if
  contains
    function check(ptr) result(q)
      ! Check if a material layer should be skipped.
      class(material), pointer :: ptr
      logical                  :: q
      if (.not. associated(ptr)) then
        q = .false.
      else
        q = ptr % lock
      end if
    end function
    subroutine top(ptr)
      ! Find the first unlocked material in the stack.
      class(material), pointer :: ptr
      ptr => this % a
      do while (check(ptr))
        ptr => ptr % material_b
      end do
    end subroutine
    subroutine next(ptr)
      ! Find the next unlocked material in the stack.
      class(material), pointer :: ptr
      if (associated(ptr)) then
        ptr => ptr % material_b
        do while (check(ptr))
          ptr => ptr % material_b
        end do
      end if
    end subroutine
    subroutine prev(ptr)
      ! Find the previous unlocked material in the stack.
      class(material), pointer :: ptr
      if (associated(ptr)) then
        ptr => ptr % material_a
        do while (check(ptr))
          ptr => ptr % material_a
        end do
      end if
    end subroutine
    subroutine bottom(ptr)
      ! Find the last unlocked material in the stack.
      class(material), pointer :: ptr
      ptr => this % b
      do while (check(ptr))
        ptr => ptr % material_a
      end do
    end subroutine
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
 
  impure function structure_count(this) result(num)
    !! Checks the number of unlocked materials in the multilayer stack.
    class(structure), target :: this
    integer                  :: num

    ! Initialize variables
    num = 0

    ! Count the number of materials
    call this % map(count)
  contains
    subroutine count
      num = num + 1
    end subroutine
  end function

  impure function structure_difference(this) result(difference)
    !! Checks how much the multilayer stack has changed recently.
    class(structure), target  :: this
    integer                   :: n, s
    real(wp)                  :: difference

    ! Initialize variables
    n = 0
    s = 0
    difference = 0

    ! Traverse all materials
    call this % map(check)

    ! Check whether we can assume convergence
    if (s == 0 .and. n <= 1) then
      difference = 0
    end if
  contains
    subroutine check(m)
      class(material) :: m
      ! Count the material types
      select type (m)
        class is (superconductor)
          s = s + 1
        class default
          n = n + 1
      end select
      ! Accumulate the difference
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
    !! Writes the density of states as a function of position and energy to a given output file.
    class(structure), target  :: this
    character(*)              :: file
    integer                   :: unit
    integer                   :: iostat
    real(wp)                  :: a, b

    ! Open output file
    open(newunit = unit, file = file, iostat = iostat, action = 'write', status = 'replace')
    if (iostat /= 0) then
      call error('Failed to open output file "' // file // '"!')
    end if

    ! Initialize variables
    b = 0

    ! Traverse all materials
    call this % map(write_density)

    ! Close output file
    close(unit = unit)
  contains
    subroutine write_density(ptr)
      class(material) :: ptr
      real(wp)        :: x
      integer         :: n, m

      ! Calculate the endpoints of this layer
      a = b
      b = b + 1/sqrt(ptr % thouless)

      ! Write out the density of states in this layer
      do m=1,size(ptr % location)
        x = a + (b-a) * ptr % location(m)
        do n=size(ptr % energy),1,-1
          ! Negative energies
          write(unit,*) x, -ptr % energy(n), ptr % density(n,m)
        end do
        do n=1,size(ptr % energy),+1
          ! Positive energies
          write(unit,*) x, +ptr % energy(n), ptr % density(n,m)
        end do
      end do
    end subroutine
  end subroutine

  impure subroutine structure_write_current(this, file)
    !! Writes the charge and spin currents as a function of position to a given output file.
    class(structure), target  :: this
    character(*)              :: file
    integer                   :: unit
    integer                   :: iostat
    real(wp)                  :: a, b

    ! Open output file
    open(newunit = unit, file = file, iostat = iostat, action = 'write', status = 'replace')
    if (iostat /= 0) then
      call error('Failed to open output file "' // file // '"!')
    end if

    ! Initialize variables
    b = 0

    ! Traverse all materials
    call this % map(write_current)

    ! Close output file
    close(unit = unit)
  contains
    subroutine write_current(ptr)
      class(material) :: ptr
      real(wp)        :: x
      integer         :: m

      ! Calculate the endpoints of this layer
      a = b
      b = b + 1/sqrt(ptr % thouless)

      ! Write out the currents in this layer
      do m=1,size(ptr % location)
        x = a + (b-a) * ptr % location(m)
        write(unit,*) x, ptr % current(:,m)
      end do
    end subroutine
  end subroutine

  impure subroutine structure_write_gap(this, file)
    !! Writes the superconducting gap as a function of position to a given output file.
    class(structure), target  :: this
    character(*)              :: file
    integer                   :: unit
    integer                   :: iostat
    real(wp)                  :: a, b

    ! Open output file
    open(newunit = unit, file = file, iostat = iostat, action = 'write', status = 'replace')
    if (iostat /= 0) then
      call error('Failed to open output file "' // file // '"!')
    end if

    ! Initialize variables
    b = 0

    ! Traverse all materials
    call this % map(write_gap)

    ! Close output file
    close(unit = unit)
  contains
    subroutine write_gap(ptr)
      class(material) :: ptr
      real(wp)        :: x
      integer         :: m

      ! Calculate the endpoints of this layer
      a = b
      b = b + 1/sqrt(ptr % thouless)

      ! Write out the gap in this layer
      do m=1,size(ptr % location)
        x = a + (b-a) * ptr % location(m)
        select type (ptr)
          class is (superconductor)
            write(unit,*) x, abs(ptr%gap(m)), atan2(im(ptr%gap(m)),re(ptr%gap(m)))/pi
          class default
            write(unit,*) x, 0.0_wp, 0.0_wp
        end select
      end do
    end subroutine
  end subroutine

  impure function structure_construct(file) result(this)
    !! Constructs a multilayer stack from a configuration file.
    type(structure)          :: this
    character(*), intent(in) :: file
    integer                  :: unit
    integer                  :: iostat
    character(len=132)       :: str, arg
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
      str = trim(adjustl(str))

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

      ! Construct a material
      i = scan(str, '[')
      j = scan(str, ']')
      if (i == 1 .and. j > i) then
        call this % push(trim(adjustl(str(i+1:j-1))))
        cycle
      end if

      ! Construct a material
      i = scan(str, ':')
      if ( i > 0 ) then
        call this % conf(trim(adjustl(str(:i-1))), trim(adjustl(str(i+1:))))
        cycle
      end if

      ! Check if this line contains garbage
      if ( str /= '' ) then
        write(str,'(i0)') line
        call error('failed to parse line ' // trim(str) // ' in the config file.')
      end if
    end do

    ! Close the config file
    close(unit = unit)

    ! Confirm that there is at least one material layer
    if (this % count() < 1) then
      call error('The material stack described by "' // file // '" has no unlocked layers!')
    end if
  end function

  !--------------------------------------------------------------------------------!
  !        @TODO: PROCEDURES FOR CALCULATING MULTILAYER INTERACTIONS               !
  !--------------------------------------------------------------------------------!

  function differential_conductance(material_a, material_b, voltage, temperature) result(conductance)
    ! Numerically calculates the differential conductance at a tunneling interface.
    class(material), intent(in) :: material_a
    class(material), intent(in) :: material_b
    real(wp),        intent(in) :: temperature
    real(wp),        intent(in) :: voltage(:)
    real(wp),       allocatable :: conductance(:)
    real(wp),       allocatable :: current(:)
    real(wp),       allocatable :: energy(:)
    real(wp),       allocatable :: dos_a(:), idos_a(:)
    real(wp),       allocatable :: dos_b(:), idos_b(:)
    integer                     :: n

    ! Allocate memory
    allocate(current(size(voltage)))
    allocate(energy(2*size(voltage)))

    allocate(dos_a(size(material_a%energy)))
    allocate(dos_b(size(material_b%energy)))
    allocate(idos_a(size(energy)))
    allocate(idos_b(size(energy)))

    ! Initialize the energy array
    call linspace(energy, 1e-6_wp, 1.50_wp)

    ! Calculate the density of states at the interface
    do n = 1,size(material_a%energy)
      dos_a(n) = material_a % propagator(n,ubound(material_a%location,1)) % get_dos()
      dos_b(n) = material_b % propagator(n,lbound(material_b%location,1)) % get_dos()
    end do

    ! Interpolate the density of states in the left material
    idos_a = interpolate(material_a%energy, dos_a, energy)

    ! Calculate the current as a function of voltage by numerical integration
    do n = 1,size(voltage)
      ! Interpolate the voltage-shifted density of states in the right material
      idos_b = interpolate(material_b%energy, dos_b, abs(energy-voltage(n)))

      ! Calculate the current for this voltage
      current(n) = integrate(energy, idos_a*idos_b*(fermi(energy-voltage(n))-fermi(energy)))
    end do

    ! Calculate the differential conductance by differentiation
    conductance = differentiate(voltage, current)

   ! Deallocate workspace memory
   deallocate(current)
   deallocate(energy)
   deallocate(dos_a)
   deallocate(dos_b)
   deallocate(idos_a)
   deallocate(idos_b)
  contains
    pure elemental function fermi(energy)
      ! Evaluates the Fermi function at the given energy and current temperature.
      real(wp), intent(in) :: energy
      real(wp)             :: fermi

      fermi = 1/(exp(energy/(temperature+1e-16))+1)
    end function
  end function
end module
