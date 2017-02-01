!> Author:   Jabir Ali Ouassou
!> Date:     2016-03-08
!> Category: Materials
!>
!> This module defines a data type 'structure', which is useful for constructing and using multilayer hybrid structures. It
!> also exports the type definitions and constructors for all class(material) types, although these should rarely be needed.

module structure_m
  use :: math_m
  use :: stdio_m
  use :: material_m
  use :: conductor_m
  use :: superconductor_m
  use :: ferromagnet_m
  use :: halfmetal_m
  private

  ! Export class(material) types
  public conductor, superconductor, ferromagnet, halfmetal

  ! Type declaration
  type, public :: structure
    class(material), pointer :: a => null()
    class(material), pointer :: b => null()
  contains
    procedure :: push                => structure_push
    procedure :: conf                => structure_conf
    procedure :: map                 => structure_map
    procedure :: gap                 => structure_gap
    procedure :: init                => structure_init
    procedure :: save                => structure_save
    procedure :: load                => structure_load
    procedure :: update              => structure_update
    procedure :: materials           => structure_materials
    procedure :: superconductors     => structure_superconductors
    procedure :: difference          => structure_difference
    procedure :: chargeviolation     => structure_chargeviolation
    procedure :: temperature         => structure_temperature
    procedure :: selfconsistency     => structure_selfconsistency
    procedure :: converge            => structure_converge
    procedure :: write_density       => structure_write_density
    procedure :: write_current       => structure_write_current
    procedure :: write_magnetization => structure_write_magnetization
    procedure :: write_gap           => structure_write_gap
  end type

  ! Type constructor
  interface structure
    module procedure structure_construct
  end interface

  ! Interface for external routines that can be mapped onto class(material) objects
  abstract interface
    subroutine mappable(ptr)
      use :: material_m
      class(material), pointer, intent(in) :: ptr
    end subroutine
  end interface

  ! Interface for external routines that can be used by convergence() calls
  abstract interface
    subroutine hook()
    end subroutine
  end interface
contains
  impure subroutine structure_push(this, string)
    !! Constructs a new class(material) object at the bottom of the multilayer stack.
    class(structure), intent(inout) :: this
    character(*),     intent(in)    :: string

    ! Construct the material
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

    ! Status information
    write(stdout,*)
    write(stdout,*) '[' // string // ']'
  contains
    impure subroutine material_allocate(ptr, str)
      !! Allocates memory for a new material layer.
      class(material), pointer, intent(inout) :: ptr
      character(*),             intent(in)    :: str

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
    character(*),     intent(in)    :: key
    character(*),     intent(in)    :: val
    character(24)                   :: str

    ! Configure the material
    if (associated(this % b)) then
      call this % b % conf(key, val)
    else
      call error('Attempted to configure a non-existant material!')
    end if

    ! Status information
    write(str,    *) key // ':'
    write(stdout, *) str // val
  end subroutine

  impure subroutine structure_map(this, routine)
    !! Maps a subroutine onto each element of the multilayer stack.
    class(structure), target  :: this
    class(material),  pointer :: ptr
    procedure(mappable)       :: routine

    ! Traverse the structure from top to bottom
    call top(ptr)
    do while (associated(ptr))
      call routine(ptr)
      call next(ptr)
    end do
  contains
    function check(ptr) result(skip)
      ! Check if a material layer should be skipped.
      class(material), pointer :: ptr
      logical                  :: skip
      if (.not. associated(ptr)) then
        skip = .false.
      else
        skip = ptr % order <= 0
      end if
    end function
    subroutine top(ptr)
      ! Find the first enabled material in the stack.
      class(material), pointer :: ptr
      ptr => this % a
      do while (check(ptr))
        ptr => ptr % material_b
      end do
    end subroutine
    subroutine next(ptr)
      ! Find the next enabled material in the stack.
      class(material), pointer :: ptr
      if (associated(ptr)) then
        ptr => ptr % material_b
        do while (check(ptr))
          ptr => ptr % material_b
        end do
      end if
    end subroutine
  end subroutine

  impure function structure_gap(this) result(gap)
    !! Obtains the mean gap in the enabled superconductor. If there are multiple such
    !! superconductors in the junction, then it returns the minimum of the mean gaps.
    class(structure), target :: this
    real(wp)                 :: gap
    logical                  :: found

    ! Initialize variables
    gap   = huge(real(wp))
    found = .false.

    ! Check the gaps
    call this % map(find)

    ! If no superconductor was found, raise an error
    if (.not. found) then
      call error('No superconductors with order > 0 in the junction!')
    end if
  contains
    subroutine find(m)
      class(material), pointer, intent(in) :: m
      select type(m)
        class is (superconductor)
          gap = min(gap, sum(abs(m%gap_function))/max(1,size(m%gap_function)))
          found = .true.
      end select
    end subroutine
  end function

  impure subroutine structure_init(this, gap)
    !! Initializes the state of the entire multilayer stack.
    class(structure), target  :: this
    complex(wp)               :: gap

    ! Initialize all material states
    call this % map(init)
  contains
    subroutine init(m)
      class(material), pointer, intent(in) :: m
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
      class(material), pointer, intent(in) :: m
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
      class(material), pointer, intent(in) :: m
      call m % load
    end subroutine
  end subroutine

  impure subroutine structure_update(this, bootstrap)
    !! Updates the state of the entire multilayer stack.
    class(structure), target   :: this
    logical,          optional :: bootstrap
    integer                    :: order

    ! Update materials in order
    do order = 1,16
      call this % map(update)
    end do
  contains
    subroutine update(m)
      class(material), pointer, intent(in) :: m

      if (m % order == order) then
        call m % update(bootstrap)
      end if
    end subroutine 
  end subroutine

  impure subroutine structure_converge(this, threshold, iterations, bootstrap, prehook, posthook)
    !! Performs a convergence procedure, where the state of every material in the stack
    !! is repeatedly updated until the residuals drop below some specified threshold 
    !! and/or a certain number of iterations have been performed. If bootstrap is set
    !! to true, the selfconsistency equations will only be solved once at the end, but
    !! not inbetween the individual iterations. If a prehook and/or posthook is given,
    !! those subroutines will be executed before/after each iteration of the update.
    class(structure), target   :: this
    real(wp),         optional :: threshold
    real(wp)                   :: threshold_
    integer,          optional :: iterations
    integer                    :: iterations_
    logical,          optional :: bootstrap
    logical                    :: bootstrap_
    procedure(hook),  optional :: prehook
    procedure(hook),  optional :: posthook
    integer                    :: materials
    integer                    :: superconductors
    integer                    :: n

    ! Set default arguments
    threshold_  = 1
    iterations_ = 0
    bootstrap_  = .false.

    ! Check optional arguments
    if (present(threshold))  threshold_  = threshold
    if (present(iterations)) iterations_ = iterations
    if (present(bootstrap))  bootstrap_  = bootstrap

    ! Count the number of materials
    materials       = this % materials()
    superconductors = this % superconductors()

    ! If we're not bootstrapping, then we have to solve the diffusion equation until convergence.
    ! If we're bootstrapping, it's only required if we have at least one enabled superconductor.
    if ((.not. bootstrap_) .or. (superconductors > 0)) then
      n = 0
      do
        ! Update counter
        n = n+1

        ! Status information
        if (bootstrap_) then
          call status_head('BOOTSTRAPPING')
        else
          call status_head('CONVERGING')
        end if
        call status_body('State difference', this % difference())
        if (.not. bootstrap_) then
          call status_body('Charge violation', this % chargeviolation())
        end if
        if (present(prehook)) then
          call prehook
        end if
        call status_body('Iteration', n)
        call status_foot

        ! Update the material state (non-selfconsistently)
        call this % update(bootstrap = bootstrap_)

        ! Write the results to files
        if (present(posthook)) then
          call posthook
        end if

        ! Exit criterion #1: one iteration is sufficient for convergence
        if ((materials == 1) .and. (superconductors == 0 .or. bootstrap_)) then
          exit
        end if

        ! Exit criterion #2: minimum number of iterations reached,
        ! and the materials converged within specified parameters
        if ((n >= iterations_) .and. (this % difference() < threshold_)) then
          exit
        end if
      end do
    end if

    ! Solve the selfconsistency equations at least once
    if (bootstrap_) then
      call this % load()
    end if
  end subroutine
 
  impure function structure_materials(this) result(num)
    !! Checks the number of enabled materials in the multilayer stack.
    class(structure), target :: this
    integer                  :: num

    ! Initialize variables
    num = 0

    ! Count the number of materials
    call this % map(count)
  contains
    subroutine count(ptr)
      class(material), pointer, intent(in) :: ptr
      num = num + 1
    end subroutine
  end function

  impure function structure_superconductors(this) result(num)
    !! Checks the number of enabled superconductors in the multilayer stack.
    class(structure), target :: this
    integer                  :: num

    ! Initialize variables
    num = 0

    ! Count the number of superconductors
    call this % map(count)
  contains
    subroutine count(ptr)
      class(material), pointer, intent(in) :: ptr
      select type (ptr)
        class is (superconductor)
          num = num + 1
      end select
    end subroutine
  end function


  impure function structure_difference(this) result(difference)
    !! Checks how much the multilayer stack has changed recently.
    class(structure), target  :: this
    real(wp)                  :: difference

    ! Traverse all materials
    call this % map(check)
  contains
    subroutine check(m)
      class(material), pointer, intent(in) :: m
      ! Accumulate the difference
      difference = max(difference, m % difference)
    end subroutine
  end function

  impure function structure_chargeviolation(this) result(difference)
    !! Checks how much the charge current varies with position. Since charge current
    !! is supposed to be conserved through the junction, this provides a measure of
    !! charge conservation violation, i.e. if the solution is physically realistic.
    class(structure), target  :: this
    real(wp)                  :: difference
    real(wp)                  :: minimum
    real(wp)                  :: maximum

    ! Set starting values
    minimum = +inf
    maximum = -inf

    ! Traverse all materials to find the most extreme currents
    call this % map(check)

    ! Calculate the difference between these extreme values
    difference = maximum - minimum
  contains
    subroutine check(m)
      class(material), pointer, intent(in) :: m
      ! Determine the charge current extrema
      maximum = max(maximum, maxval(m % current(0,:)))
      minimum = min(minimum, minval(m % current(0,:)))
    end subroutine
  end function

  impure subroutine structure_selfconsistency(this, scheme)
    !! Controls the selfconsistency scheme used for superconductors in the stack.
    class(structure), target  :: this
    integer                   :: scheme

    call this % map(set)
  contains
    subroutine set(m)
      class(material), pointer, intent(in) :: m
      select type (m)
        class is (superconductor)
          m % selfconsistency = scheme
      end select
    end subroutine
  end subroutine

  impure subroutine structure_temperature(this, temperature)
    !! Modifies the temperature of the multilayer stack.
    class(structure), target  :: this
    real(wp)                  :: temperature

    ! Update the temperature of all materials
    call this % map(set)
  contains
    subroutine set(m)
      class(material), pointer, intent(in) :: m
      m % temperature = temperature
    end subroutine
  end subroutine

  impure subroutine structure_write_density(this, file)
    !! Writes the density of states as a function of position and energy to a given output file.
    class(structure), target  :: this
    character(*)              :: file
    integer                   :: unit
    real(wp)                  :: a, b

    ! Open output file
    unit = output(file)

    ! Initialize variables
    b = 0

    ! Write out the header line
    write(unit,'(*(a,:,"	"))') '# Position          ', &
                                '  Energy            ', &
                                '  Density of states '

    ! Traverse all materials
    call this % map(write_density)

    ! Close output file
    close(unit = unit)
  contains
    subroutine write_density(ptr)
      class(material), pointer, intent(in) :: ptr
      real(wp)                             :: x
      integer                              :: n, m

      ! Calculate the endpoints of this layer
      a = b
      b = b + 1/sqrt(ptr % thouless)

      ! Write out the density of states in this layer
      do m=1,size(ptr % location)
        x = a + (b-a) * ptr % location(m)
        do n=size(ptr % energy),1,-1
          ! Negative energies
          write(unit,'(*(es20.12e3,:,"	"))') x, -ptr % energy(n), ptr % density(n,m)
        end do
        do n=1,size(ptr % energy),+1
          ! Positive energies
          write(unit,'(*(es20.12e3,:,"	"))') x, +ptr % energy(n), ptr % density(n,m)
        end do
      end do
    end subroutine
  end subroutine

  impure subroutine structure_write_current(this, file)
    !! Writes the charge and spin currents as a function of position to a given output file.
    class(structure), target  :: this
    character(*)              :: file
    integer                   :: unit
    real(wp)                  :: a, b

    ! Open output file
    unit = output(file)

    ! Initialize variables
    b = 0

    ! Write out the header line
    write(unit,'(*(a,:,"	"))') '# Position          ', &
                                '  Charge current    ', &
                                '  Spin-x current    ', &
                                '  Spin-y current    ', &
                                '  Spin-z current    '

    ! Traverse all materials
    call this % map(write_current)

    ! Close output file
    close(unit = unit)
  contains
    subroutine write_current(ptr)
      class(material), pointer, intent(in) :: ptr
      real(wp)                             :: x
      integer                              :: m

      ! Calculate the endpoints of this layer
      a = b
      b = b + 1/sqrt(ptr % thouless)

      ! Write out the currents in this layer
      do m=1,size(ptr % location)
        x = a + (b-a) * ptr % location(m)
        write(unit,'(*(es20.12e3,:,"	"))') x, ptr % current(:,m)
      end do
    end subroutine
  end subroutine

  impure subroutine structure_write_magnetization(this, file)
    !! Writes the induced magnetization as a function of position to a given output file.
    class(structure), target  :: this
    character(*)              :: file
    integer                   :: unit
    real(wp)                  :: a, b

    ! Open output file
    unit = output(file)

    ! Initialize variables
    b = 0

    ! Write out the header line
    write(unit,'(*(a,:,"	"))') '# Position          ', &
                                '  Mx magnetization  ', &
                                '  My magnetization  ', &
                                '  Mz magnetization  '

    ! Traverse all materials
    call this % map(write_magnetization)

    ! Close output file
    close(unit = unit)
  contains
    subroutine write_magnetization(ptr)
      class(material), pointer, intent(in) :: ptr
      real(wp)                             :: x
      integer                              :: m

      ! Calculate the endpoints of this layer
      a = b
      b = b + 1/sqrt(ptr % thouless)

      ! Write out the magnetization in this layer
      do m=1,size(ptr % location)
        x = a + (b-a) * ptr % location(m)
        write(unit,'(*(es20.12e3,:,"	"))') x, ptr % magnetization(:,m)
      end do
    end subroutine
  end subroutine

  impure subroutine structure_write_gap(this, file)
    !! Writes the superconducting gap as a function of position to a given output file.
    class(structure), target  :: this
    character(*)              :: file
    integer                   :: unit
    real(wp)                  :: a, b

    ! Open output file
    unit = output(file)

    ! Initialize variables
    b = 0

    ! Write out the header line
    write(unit,'(*(a,:,"	"))') '# Position          ', &
                                '  Gap magnitude     ', &
                                '  Gap phase         '
    ! Traverse all materials
    call this % map(write_gap)

    ! Close output file
    close(unit = unit)
  contains
    subroutine write_gap(ptr)
      class(material), pointer, intent(in) :: ptr
      real(wp)                             :: x
      real(wp)                             :: p
      integer                              :: m

      ! Calculate the endpoints of this layer
      a = b
      b = b + 1/sqrt(ptr % thouless)

      ! Write out the gap in this layer
      do m=1,size(ptr % location)
        p = ptr % location(m)
        x = a + (b-a)*p
        select type (ptr)
          class is (superconductor)
            write(unit,'(*(es20.12e3,:,"	"))') x, abs(ptr%get_gap(p)), atan2(im(ptr%get_gap(p)),re(ptr%get_gap(p)))/pi
          class default
            write(unit,'(*(es20.12e3,:,"	"))') x, 0.0_wp, 0.0_wp
        end select
      end do
    end subroutine
  end subroutine

  impure function structure_construct(file) result(this)
    !! Constructs a multilayer stack from a configuration file.
    use :: stdio_m

    type(structure)          :: this
    character(*), intent(in) :: file
    integer                  :: unit
    integer                  :: iostat
    character(len=2048)      :: str, arg
    integer                  :: line, i, j, k

    ! Initialize variables
    line   = 0
    unit   = 0
    iostat = 0

    ! Open the config file
    unit = input(file)

    ! Command arguments
    call get_command(str)
    write(stdout, '(a)') trim(str)

    ! Status information
    call status_box('CONFIGURATION')

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
          call error('Missing command line arguments.')
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
        call error('Failed to parse line ' // trim(str) // ' in the config file.')
      end if
    end do

    ! Close the config file
    close(unit = unit)

    ! Confirm that there is at least one material layer
    if (this % materials() < 1) then
      call error('The material stack described by "' // file // '" has no layers with order > 0!')
    end if
  end function
end module
