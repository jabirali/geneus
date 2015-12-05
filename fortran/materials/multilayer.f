! This module...
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
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-12-05
! Updated: 2015-12-05

module mod_multilayer
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
    string = "-"
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

    ! Configure the material layers
    s = 0
    f = 0
    c = 0
    this => m % top
    do while (associated(this))
      select type(this)
        class is (superconductor)
          s = s + 1

          m % s(s) % magnetization_a = [0,0,0]
          m % s(s) % magnetization_b = [0,0,0]

          write(*,'(a)') '[interface]'
          call config(unit, 'superconductor', 'transparent',   s, -1, m % s(s) % transparent_a)
          call config(unit, 'superconductor', 'reflecting',    s, -1, m % s(s) % reflecting_a)
          call config(unit, 'superconductor', 'conductance',   s, -1, m % s(s) % conductance_a)
          call config(unit, 'superconductor', 'spinmixing',    s, -1, m % s(s) % spinmixing_a)
          call config(unit, 'superconductor', 'polarization',  s, -1, m % s(s) % polarization_a)
          call config(unit, 'superconductor', 'magnetization', s, -1, m % s(s) % magnetization_a)
          write(*,*)

          write(*,'(a)') '[superconductor]'
          write(*,*)

          write(*,'(a)') '[interface]'
          call config(unit, 'superconductor', 'transparent',   s, +1, m % s(s) % transparent_b)
          call config(unit, 'superconductor', 'reflecting',    s, +1, m % s(s) % reflecting_b)
          call config(unit, 'superconductor', 'conductance',   s, +1, m % s(s) % conductance_b)
          call config(unit, 'superconductor', 'spinmixing',    s, +1, m % s(s) % spinmixing_b)
          call config(unit, 'superconductor', 'polarization',  s, +1, m % s(s) % polarization_b)
          call config(unit, 'superconductor', 'magnetization', s, +1, m % s(s) % magnetization_b)
          write(*,*)

        class is (ferromagnet)
          f = f + 1

          m % f(f) % magnetization_a = [0,0,0]
          m % f(f) % magnetization_b = [0,0,0]

          write(*,'(a)') '[interface]'
          call config(unit, 'ferromagnet',    'transparent',   f, -1, m % f(f) % transparent_a)
          call config(unit, 'ferromagnet',    'reflecting',    f, -1, m % f(f) % reflecting_a)
          call config(unit, 'ferromagnet',    'conductance',   f, -1, m % f(f) % conductance_a)
          call config(unit, 'ferromagnet',    'spinmixing',    f, -1, m % f(f) % spinmixing_a)
          call config(unit, 'ferromagnet',    'polarization',  f, -1, m % f(f) % polarization_a)
          call config(unit, 'ferromagnet',    'magnetization', f, -1, m % f(f) % magnetization_a)
          write(*,*)

          write(*,'(a)') '[ferromagnet]'
          write(*,*)

          write(*,'(a)') '[interface]'
          call config(unit, 'ferromagnet',    'transparent',   f, +1, m % f(f) % transparent_b)
          call config(unit, 'ferromagnet',    'reflecting',    f, +1, m % f(f) % reflecting_b)
          call config(unit, 'ferromagnet',    'conductance',   f, +1, m % f(f) % conductance_b)
          call config(unit, 'ferromagnet',    'spinmixing',    f, +1, m % f(f) % spinmixing_b)
          call config(unit, 'ferromagnet',    'polarization',  f, +1, m % f(f) % polarization_b)
          call config(unit, 'ferromagnet',    'magnetization', f, +1, m % f(f) % magnetization_b)
          write(*,*)

        class is (conductor)
          c = c + 1

          m % c(c) % magnetization_a = [0,0,0]
          m % c(c) % magnetization_b = [0,0,0]

          write(*,'(a)') '[interface]'
          call config(unit, 'conductor',      'transparent',   c, -1, m % c(c) % transparent_a)
          call config(unit, 'conductor',      'reflecting',    c, -1, m % c(c) % reflecting_a)
          call config(unit, 'conductor',      'conductance',   c, -1, m % c(c) % conductance_a)
          call config(unit, 'conductor',      'spinmixing',    c, -1, m % c(c) % spinmixing_a)
          call config(unit, 'conductor',      'polarization',  c, -1, m % c(c) % polarization_a)
          call config(unit, 'conductor',      'magnetization', c, -1, m % c(c) % magnetization_a)
          write(*,*)

          write(*,'(a)') '[conductor]'
          write(*,*)

          write(*,'(a)') '[interface]'
          call config(unit, 'conductor',      'transparent',   c, +1, m % c(c) % transparent_b)
          call config(unit, 'conductor',      'reflecting',    c, +1, m % c(c) % reflecting_b)
          call config(unit, 'conductor',      'conductance',   c, +1, m % c(c) % conductance_b)
          call config(unit, 'conductor',      'spinmixing',    c, +1, m % c(c) % spinmixing_b)
          call config(unit, 'conductor',      'polarization',  c, +1, m % c(c) % polarization_b)
          call config(unit, 'conductor',      'magnetization', c, +1, m % c(c) % magnetization_b)
          write(*,*)
      end select

      this => this % material_b
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
