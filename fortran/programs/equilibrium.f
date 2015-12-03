! This program calculates the equilibrium propagators in a superconducting multilayer structure, and saves
! some physical observables such as the density of states and superconducting gap to separate output files.
!
! WORK IN PROGRESS:
!   The development of this program was prompted by the maintenance issues with density.f,
!   and this propgram is supposed to generalize and eventually replace density.f altogether.
!   The ideal is to make it run based on a config file (using mod_config) instead of command
!   line options, and to eventually reduce code duplication between density.f and critical.f.
!
!   Hopefully, the multilayer structure defined here will also replace and supersede the
!   methods in hybrid.f, and make it easier to construct and prototype new driver programs
!   based on superconducting materials. As well as making it easily extensible for new materials.
!   (This should of course be moved to a separate file when the program works.)
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-12-03
! Updated: 2015-12-03

module mod_local
  use mod_hybrid
  use mod_config
  implicit none

  type :: multilayer
    class(material),      pointer     :: top
    type(superconductor), allocatable :: s(:)
    type(ferromagnet),    allocatable :: f(:)
    type(conductor),      allocatable :: c(:)
  contains
    procedure :: config => multilayer_config
    procedure :: update => multilayer_update
  end type
contains
  impure subroutine multilayer_config(m, unit)
    ! This subroutine initializes the multilayer based on a config file.
    class(multilayer), target, intent(inout) :: m
    integer,                   intent(in)    :: unit
    integer                                  :: iostat

    character(len=132)               :: string
    integer                          :: s
    integer                          :: f
    integer                          :: c

    real(wp), allocatable :: energies(:)
    integer :: n

    class(material), pointer :: prev
    class(material), pointer :: this

    n = 800
    call config(unit, 'global', 'energies', 1, 0, n)
    allocate(energies(n))
    call energy_range(energies)


    ! Count the material layers
    rewind(unit)
    iostat = 0
    string = ""
    s = 0
    f = 0
    c = 0
    do while (iostat == 0)
      select case (string)
        case ('[superconductor]')
          s = s + 1
        case ('[ferromagnet]')
          f = f + 1
        case ('[conductor]')
          c = c + 1
      end select
      call config(unit, iostat, string)
    end do

    ! Make sure that there is at least one layer
    if (s == 0) then
      call error('at least one superconducting layer is required!')
    end if

    ! Allocate memory for the layers
    allocate(m % s(s))
    allocate(m % f(f))
    allocate(m % c(c))

    ! Initialize the material layers
    rewind(unit)
    iostat = 0
    string = ""
    prev => null()
    this => null()
    s = 0
    f = 0
    c = 0
    do while (iostat == 0)
      select case (string)
        case ('[superconductor]')
          s = s + 1
          prev => this
          this => m % s(s)

          m % s(s) = superconductor(energies)
        case ('[ferromagnet]')
          f = f + 1
          prev => this
          this => m % f(f)

          m % f(f) = ferromagnet(energies)
        case ('[conductor]')
          c = c + 1
          prev => this
          this => m % c(c)

          m % c(c) = conductor(energies)
      end select
      if (associated(this)) then
        if (associated(prev)) then
          prev % material_b => this
          this % material_a => prev
        else
          m % top => this
        end if
      end if
      call config(unit, iostat, string)
    end do
  end subroutine

  impure subroutine multilayer_update(m)
    ! This subroutine updates the state of the entire hybrid structure.
    class(multilayer), intent(in) :: m
    class(material),   pointer    :: p
    integer                       :: n

    ! Initialize the material pointer to the top of the stack
    p => m % top

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

program equilibrium
  use mod_stdio,  only: error
  use mod_config, only: config
  use mod_hybrid
  use mod_local
  implicit none

  ! Setting that control the program
  integer :: unit            = 0
  integer :: iostat          = 0
  integer :: energies        = 800
  integer :: debug           = 0

  ! Computational model of the system
  real(wp), allocatable :: e(:)
  type(multilayer), target :: system

  ! Open configuration file
  open(newunit=unit, file='config.ini', action='read', status='old', iostat=iostat)
  if (iostat /= 0) then
    call error('failed to open configuration file ''config.ini''!')
  end if

  ! Process configuration file
  write(*,'(a)') '[equilibrium]'
  call config(unit, 'equilibrium', 'debug',    1, 0, debug)

  ! Count the material layers

  ! Construct material layers

  ! Computation
  call system % config(unit)
  call system % update
  call system % update

  ! Close configuration file
  close(unit)
end program
