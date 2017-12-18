!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This module defines a data type 'structure', which is useful for constructing and using multilayer hybrid structures. It
!> also exports the type definitions and constructors for all class(material) types, although these should rarely be needed.

module structure_m
  use :: stdio_m
  use :: condmat_m
  use :: material_m
  use :: conductor_m
  use :: superconductor_m
  use :: ferromagnet_m
  use :: halfmetal_m
  private

  ! Export class(material) types
  public material, conductor, superconductor, ferromagnet, halfmetal

  ! Type declaration
  type, public :: structure
    ! Endpoints of the contained linked list
    class(material), pointer :: a => null()   !! First material
    class(material), pointer :: b => null()   !! Last  material

    ! Output units for physical observables
    integer, allocatable :: supercurrent      !! Output unit (allocate to write supercurrents to file)
    integer, allocatable :: lossycurrent      !! Output unit (allocate to write lossycurrents to file)
    integer, allocatable :: accumulation      !! Output unit (allocate to write accumulations to file)
    integer, allocatable :: correlation       !! Output unit (allocate to write correlations  to file)
    integer, allocatable :: density           !! Output unit (allocate to write density of states to file)
  contains
    ! Basic construction and management
    procedure :: push                => structure_push                !! Construct a single layer
    procedure :: conf                => structure_conf                !! Configure a single layer
    procedure :: cmap                => structure_cmap                !! Configure  all layers
    procedure :: fmap                => structure_fmap                !! Manipulate all layers

    ! Manipulation of the physical state
    procedure :: initialize          => structure_initialize          !! Reset the physical state
    procedure :: save                => structure_save                !! Save the physical state
    procedure :: load                => structure_load                !! Load the physical state
    procedure :: update              => structure_update              !! Update the physical state
    procedure :: update_prehook      => structure_update_prehook      !! Execute all update prehooks
    procedure :: update_posthook     => structure_update_posthook     !! Execute all update posthooks
    procedure :: converge            => structure_converge            !! Update until convergence
    procedure :: write               => structure_write               !! Write out observables

    ! Auxiliary helper functions
    procedure :: difference          => structure_difference          !! Check how much the physical state changes
    procedure :: materials           => structure_materials           !! Check the number of enabled materials
    procedure :: superconductors     => structure_superconductors     !! Check the number of enables superconductors
    procedure :: chargeviolation     => structure_chargeviolation     !! Check the violation of charge conservation
    procedure :: gap                 => structure_gap                 !! Check the minimum superconducting gap
  end type

  ! Type constructor
  interface structure
    module procedure structure_construct
  end interface

  ! Interface for external routines that can be mapped onto class(material) objects
  abstract interface
    subroutine mappable(ptr)
      use :: material_m
      class(material), pointer :: ptr
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
      call alloc(this % b, string)
      call this % b % construct
      this % a => this % b
    else
      ! This is not the first layer in the structure
      call alloc(this % b % material_b, string)
      call this % b % material_b % construct
      this % b % material_b % material_a => this % b
      this % b => this % b % material_b
    end if

    ! Status information
    write(stdout,*)
    write(stdout,*) '[' // string // ']'
  contains
    subroutine alloc(ptr, str)
      !! Allocates memory for a new material layer.
      class(material), pointer :: ptr
      character(*)             :: str

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
  end subroutine

  impure subroutine structure_conf(this, key, val)
    !! Configures the last material pushed to the multilayer stack.
    !!
    !! @TODO
    !!   Add global config options for the entire stack.
    class(structure), intent(inout) :: this
    character(*),     intent(in)    :: key
    character(*),     intent(in)    :: val
    character(24)                   :: str

    ! Status information
    write(str,    *) key // ':'
    write(stdout, *) str // val

    ! Configuration procedure
    if (.not. associated(this % b)) then
      ! Global configuration since no materials exist
      select case(key)
        case default
          call warning("Unknown structure option '" // key // "' ignored.")
      end select
    else
      ! Local configuration of only the last material
      call this % b % conf(key, val)
    end if
  end subroutine

  impure subroutine structure_cmap(this, key, val)
    !! Maps a configuration option onto each element of the multilayer stack.
    class(structure), intent(inout) :: this
    character(*),     intent(in)    :: key
    class(*),         intent(in)    :: val
    character(1024)                 :: str

    ! Parse the config option
    str = ''
    select type(val)
      type is (character(*))
        str = val
      type is (logical)
        write(str, '(l1)') val
      type is (integer)
        write(str, '(i0)') val
      type is (real(wp))
        write(str, '(g0)') val
    end select

    ! Map it onto each material
    call this % fmap(conf)
  contains
    subroutine conf(m)
      class(material), pointer :: m

      ! Configure the material
      call m % conf(key, trim(adjustl(str)))
    end subroutine
  end subroutine

  impure subroutine structure_fmap(this, routine)
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

      ! Decide whether to skip this layer
      if (.not. associated(ptr)) then
        ! There are no layers left to skip
        skip = .false.
      else
        ! Decide based on the user config
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
    call this % fmap(find)

    ! If no superconductor was found, raise an error
    if (.not. found) then
      call error('No superconductors with order > 0 in the junction!')
    end if
  contains
    subroutine find(m)
      class(material), pointer :: m
      select type(m)
        class is (superconductor)
          gap = min(gap, sum(abs(m%gap_function))/max(1,size(m%gap_function)))
          found = .true.
      end select
    end subroutine
  end function

  impure subroutine structure_initialize(this)
    !! Initializes the state of the entire multilayer stack.
    class(structure), target  :: this

    ! Initialize all material states
    call this % fmap(initialize)
  contains
    subroutine initialize(m)
      class(material), pointer :: m
      call m % initialize
    end subroutine
  end subroutine

  impure subroutine structure_save(this)
    !! Saves the state of the entire multilayer stack.
    class(structure), target  :: this

    ! Save all material states
    call this % fmap(save)
  contains
    subroutine save(m)
      class(material), pointer :: m
      call m % save
    end subroutine
  end subroutine

  impure subroutine structure_load(this)
    !! Loads the saved state of the multilayer stack.
    class(structure), target  :: this

    ! Load all material states
    call this % fmap(load)
  contains
    subroutine load(m)
      class(material), pointer :: m
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
      call this % fmap(update)
    end do
  contains
    subroutine update(m)
      class(material), pointer :: m

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

    ! If we're not bootstrapping, then we wish to execute all posthook actions at least once
    if (.not. bootstrap_) then
      call this % update_posthook
    end if

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
        call this % write

        ! Extra actions defined by the user
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
  end subroutine

  impure subroutine structure_update_prehook(this)
    !! Execute all update prehooks.
    class(structure) :: this

    ! Traverse all materials
    call this % fmap(prehook)
  contains
    subroutine prehook(ptr)
      class(material), pointer :: ptr

      call ptr % update_prehook
    end subroutine
  end subroutine

  impure subroutine structure_update_posthook(this)
    !! Execute all update posthooks.
    class(structure) :: this

    ! Traverse all materials
    call this % fmap(posthook)
  contains
    subroutine posthook(ptr)
      class(material), pointer :: ptr

      call ptr % update_posthook
    end subroutine
  end subroutine
 
  impure function structure_materials(this) result(num)
    !! Checks the number of enabled materials in the multilayer stack.
    class(structure), target :: this
    integer                  :: num

    ! Initialize variables
    num = 0

    ! Count the number of materials
    call this % fmap(count)
  contains
    subroutine count(ptr)
      class(material), pointer :: ptr
      num = num + 1
    end subroutine
  end function

  impure function structure_superconductors(this) result(num)
    !! Checks the number of selfconsistent superconductors in the multilayer stack.
    class(structure), target :: this
    integer                  :: num

    ! Initialize variables
    num = 0

    ! Count the number of superconductors
    call this % fmap(count)
  contains
    subroutine count(ptr)
      class(material), pointer :: ptr
      select type (ptr)
        class is (superconductor)
          if (ptr % selfconsistency > 0) then
            num = num + 1
          end if
      end select
    end subroutine
  end function

  impure function structure_difference(this) result(difference)
    !! Checks how much the multilayer stack has changed recently.
    class(structure), target  :: this
    real(wp)                  :: difference

    ! Initialization
    difference = 0

    ! Traverse all materials
    call this % fmap(check)
  contains
    subroutine check(m)
      class(material), pointer :: m
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
    call this % fmap(check)

    ! Calculate the difference between these extreme values
    difference = maximum - minimum
  contains
    subroutine check(m)
      class(material), pointer :: m
      ! Determine the charge current extrema
      maximum = max(maximum, maxval(m % supercurrent(0,:) + m % lossycurrent(0,:)))
      minimum = min(minimum, minval(m % supercurrent(0,:) + m % lossycurrent(0,:)))
    end subroutine
  end function

  impure subroutine structure_write(this)
    !! Writes physical observables to output files.
    class(structure), target  :: this
    real(wp)                  :: a, b

    ! Initialize variables
    b = 0

    ! Traverse all materials
    call this % fmap(writer)

    ! Sync data to output files
    call sync(this % accumulation)
    call sync(this % supercurrent)
    call sync(this % lossycurrent)
    call sync(this % correlation)
    call sync(this % density)
  contains
    subroutine sync(unit)
      integer, allocatable, intent(in) :: unit

      if (allocated(unit)) then
        flush(unit)
        rewind(unit)
      end if
    end subroutine

    subroutine writer(ptr)
      class(material), pointer :: ptr
      real(wp)                 :: p, z
      integer                  :: n, m

      ! Calculate the endpoints
      a = b
      b = b + ptr % length

      ! Loop over all positions
      do m=1,size(ptr % location)
        ! Calculate current position
        p = ptr % location(m)
        z = a + (b-a)*p

        ! Write accumulations
        if (allocated(this % accumulation)) then
          write(this % accumulation, '(*(es20.12e3,:,"	"))') &
            z, ptr % accumulation(:,m)
        end if

        ! Write supercurrents
        if (allocated(this % supercurrent)) then
          write(this % supercurrent, '(*(es20.12e3,:,"	"))') &
            z, ptr % supercurrent(:,m)
        end if

        ! Write dissipative currents
        if (allocated(this % lossycurrent)) then
          write(this % lossycurrent, '(*(es20.12e3,:,"	"))') &
            z, ptr % lossycurrent(:,m)
        end if

        ! Write superconducting correlations
        if (allocated(this % correlation)) then
          write(this % correlation,'(*(es20.12e3,:,"	"))') &
            z, abs(ptr % correlation(m)), arg(ptr % correlation(m))/pi
        end if

        ! Write density of states
        if (allocated(this % density)) then
          ! Negative energies
          do n=size(ptr % energy),1,-1
            write(this % density,'(*(es20.12e3,:,"	"))') &
              z, -ptr % energy(n), ptr % density(n,m,4:7)
          end do
          ! Positive energies
          do n=1,size(ptr % energy),+1
            write(this % density,'(*(es20.12e3,:,"	"))') &
              z, +ptr % energy(n), ptr % density(n,m,0:3)
          end do
        end if
      end do
    end subroutine
  end subroutine

  impure function structure_construct() result(this)
    !! Constructs a multilayer stack from a configuration file.

    type(structure)          :: this
    character(len=4096)      :: file
    integer                  :: unit
    integer                  :: iostat
    integer                  :: status
    character(len=2048)      :: str, arg
    integer                  :: line, i, j, k

    ! Initialize variables
    line   = 0
    unit   = 0
    iostat = 0
    status = 0

    ! Open the config file
    call get_command_argument(1, file, status=status)
    if ( status /= 0 ) then
      call error('Missing the command line argument #1, i.e. the name of configuration file.')
    end if
    unit = input(trim(file))

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
        read(str(i+1:j-1), *, iostat = status) k
        if ( status /= 0 ) then
          call error('Failed to parse parameter ' // str(i:j) // ' defined by the configuration file.')
        end if
        call get_command_argument(k+1, arg, status=status)
        if ( status /= 0 ) then
          write(arg,'(i0)') k+1
          call error('Missing command line argument #' // trim(arg) // ', i.e. parameter ' // str(i:j) // ' in the config file.')
        end if
        str = str(:i-1) // trim(arg) // str(j+1:)
        i = scan(str, '{')
        j = scan(str, '}')
      end do

      ! Construct the material
      i = scan(str, '[')
      j = scan(str, ']')
      if (i == 1 .and. j > i) then
        call this % push(trim(adjustl(str(i+1:j-1))))
        cycle
      end if

      ! Configure the material
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

    ! Initialize all materials
    call this % initialize
  end function
end module
