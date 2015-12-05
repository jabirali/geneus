! This module define the interface 'config', which can be used to read and parse configuration settings from file or command line.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-12-03
! Updated: 2015-12-03

module mod_config 
  use mod_stdio, only: error
  use mod_math,  only: wp
  implicit none
  private

  interface config
    module procedure config_scalar, config_vector, config_next
  end interface

  public config
contains
  impure subroutine config_readline(unit, iostat, output, reverse)
    ! This subroutine reads one non-empty line from an initialization file, strips it of comments and whitespace, and returns it.
    ! If an error occurs during reading, such as an end-of-file or read error, this will be reflected by a non-zero iostat value.
    integer,            intent(in)  :: unit
    integer,            intent(out) :: iostat
    logical, optional,  intent(in)  :: reverse
    character(len=132), intent(out) :: output
    character(len=396)              :: input
    character(len=132)              :: string
    integer                         :: n, m, k
    integer                         :: conts

    ! Variable initialization
    m      = 1
    iostat = 0
    output = ""
    input  = ""
    conts  = 0

    ! Process one non-empty line before returning
    do while (output == "")
      ! Move backwards if 'reverse' is enabled
      if (present(reverse)) then
        if (reverse) then
          backspace(unit)
          backspace(unit)
        end if
      end if

      ! Read one line from file
      read(unit,'(a)',iostat=iostat) input
      if (iostat /= 0) then
        exit
      end if

      ! Skip lines without sections or statements
      if (scan(input,'[=]') == 0) then
        cycle
      end if

      ! Process the input line
      do n=1,len(input)
        select case(input(n:n))
          ! Whitespace will be ignored
          case (' ', char(11))
            cycle

          ! '#' is interpreted as a start of a comment
          case ('#')
            exit

          ! '$' is interpreted as a command line argument
          case ('$')
            call get_command_argument(number=1, value=string, length=k)
            if (k <= 0) then
              call error("failed to parse parts of the configuration file due to missing command line arguments.")
            end if

          ! ',' may be interpreted as a line continuation
          case (',')
            if (input(n+1:) == '') then
              conts = conts + 1
              read(unit,'(a)',iostat=iostat) input(n+1:)
              if (iostat /= 0) then
                call error("failed to parse parts of the configuration file due to a bad line continuation.")
              end if
            end if
        end select

        ! Copy characters to output string
        output(m:m) = input(n:n) 
        m = m + 1
      end do
    end do

    ! Revert continuations by moving back
    do n=1,conts
      backspace(unit)
    end do
  end subroutine

  impure subroutine config_section(unit, iostat, section)
    ! This subroutine searches an initialization file for a given section, and returns an iostat:
    ! iostat = 0 if section found, iostat < 0 if end of file reached, iostat > 0 for misc errors.
    integer,            intent(in)  :: unit
    integer,            intent(out) :: iostat
    character(len= * ), intent(in)  :: section
    character(len=132)              :: string

    ! File and variable initialization
    string = ""
    iostat = 0

    ! Loop until the section is found or an error occurs
    do while (iostat == 0)
      ! Exit if the correct section is found
      if (string == '['//section//']') then
        exit
      end if

      ! Read another line from the file
      call config_readline(unit, iostat, string)
    end do
  end subroutine

  impure subroutine config_setting(unit, iostat, setting, value)
    ! This subroutine searches an initialization section for a given setting, and returns an iostat:
    ! iostat = 0 if setting found, iostat < 0 if end of section reached, iostat > 0 for misc errors.
    ! If the setting is successfully found, its value will be returned through the argument 'value'.
    integer,            intent(in)  :: unit
    integer,            intent(out) :: iostat
    character(len= * ), intent(in)  :: setting
    character(len=132), intent(out) :: value
    character(len=132)              :: string
    integer                         :: n

    ! Variable initialization
    value  = ""
    string = ""
    iostat = 0
    n      = 0

    ! Loop until the setting is found or an error occurs
    do while (iostat == 0)
      ! If a new section is reached, abort with a negative iostat
      if (string(1:1) == '[') then
        iostat = -2
        exit
      end if

      ! If this is the setting we are looking for, return its value
      n = scan(string, '=')
      if (string(1:n-1) == setting) then
        value = string(n+1:132)
        exit
      end if

      ! Read another line from the input file
      call config_readline(unit, iostat, string)
    end do
  end subroutine

  impure subroutine config_interface(unit, iostat, region)
    ! This subroutine moves to the next [interface] section of the initialization file if region>0, or the previous one if region<0.
    integer,            intent(in)  :: unit
    integer,            intent(in)  :: region
    integer,            intent(out) :: iostat
    character(len=132)              :: string
    integer                         :: n

    iostat = 0
    string = ""
    if (region < 0) then
      ! Loop until the interface is found or an error occurs
      do while (iostat == 0)
        ! Exit if a new section is found
        if (string(1:1) == '[') then
          if (string /= '[interface]') then
            iostat = -2
          end if
          exit
        end if

        ! Read another line from the file
        call config_readline(unit, iostat, string, reverse=.true.)
      end do
    else if (region > 0) then
      ! Loop until the interface is found or an error occurs
      do while (iostat == 0)
        ! Exit if a new section is found
        if (string(1:1) == '[') then
          if (string /= '[interface]') then
            iostat = -2
          end if
          exit
        end if

        ! Read another line from the file
        call config_readline(unit, iostat, string)
      end do
    end if
  end subroutine

  impure subroutine config_scalar(unit, section, setting, number, region, variable)
    ! This subroutine searches an initialization file for a section and setting, and initializes a scalar variable.
    character(len= * ), intent(in)    :: section
    character(len= * ), intent(in)    :: setting
    class(*),           intent(inout) :: variable
    integer,            intent(in)    :: number
    integer,            intent(in)    :: region
    integer,            intent(in)    :: unit
    integer                           :: iostat
    character(len=132)                :: string

    ! Check the global section first
    call string_check('global', 1, 0)
    call string_parse

    ! Check the local section second
    if (number > 0) then
      call string_check(section, number, region)
      call string_parse
    end if

    ! Write out the current value
    call string_print

    ! Reset the input file
    rewind(unit)
  contains
    subroutine string_check(section, number, region)
      ! Check for a variable in a given section
      character(len=*) :: section
      integer          :: number
      integer          :: region
      integer          :: n

      string = ""
      iostat = 0
      do n=1,number
        call config_section(unit, iostat, section)
        if (iostat /= 0) then
          exit
        end if
      end do
      if (iostat == 0) then
        call config_interface(unit, iostat, region)
      end if
      if (iostat == 0) then
        call config_setting(unit, iostat, setting, string)
      end if
    end subroutine

    subroutine string_parse
      ! If the string is set, parse it
      if (string /= "") then
        select type(variable)
          type is (logical)
            read(string, '(l1)',    iostat=iostat) variable
          type is (integer)
            read(string, '(i10)',   iostat=iostat) variable
          type is (real(wp))
            read(string, '(g24.0)', iostat=iostat) variable
          type is (character(*))
            read(string, '(a)',     iostat=iostat) variable
        end select
        if (iostat /= 0) then
          call error("failed to parse '" // setting // '=' // trim(string) // "' as a scalar.")
        end if
      end if
    end subroutine

    subroutine string_print
      ! Print a config option
      write(string, '(a)') trim(setting)

      select type(variable)
        type is (logical)
          write(string(21:), '("= ",l10)')   variable
        type is (integer)
          write(string(21:), '("= ",i10)')   variable
        type is (real(wp))
          write(string(21:), '("= ",f10.5)') variable
      end select

      write(*, '(a)') string
    end subroutine
  end subroutine

  impure subroutine config_vector(unit, section, setting, number, region, variable)
    ! This subroutine searches an initialization file for a section and setting, and initializes a vector variable.
    character(len= * ), intent(in)    :: section
    character(len= * ), intent(in)    :: setting
    class(*),           intent(inout) :: variable(:)
    integer,            intent(in)    :: number
    integer,            intent(in)    :: region
    integer,            intent(in)    :: unit
    integer                           :: iostat
    character(len=132)                :: string

    ! Check the global section first
    call string_check('global', 1, 0)
    call string_parse

    ! Check the local section second
    if (number > 0) then
      call string_check(section, number, region)
      call string_parse
    end if

    ! Write out the current value
    call string_print

    ! Reset the input file
    rewind(unit)
  contains
    subroutine string_check(section, number, region)
      ! Check for a variable in a given section
      character(len=*) :: section
      integer          :: number
      integer          :: region
      integer          :: n

      string = ""
      iostat = 0
      do n=1,number
        call config_section(unit, iostat, section)
        if (iostat /= 0) then
          exit
        end if
      end do
      if (iostat == 0) then
        call config_interface(unit, iostat, region)
      end if
      if (iostat == 0) then
        call config_setting(unit, iostat, setting, string)
      end if
    end subroutine

    subroutine string_parse
      ! If the string is set, parse it
      if (string /= "") then
        select type(variable)
          type is (real(wp))
            read(string, *, iostat=iostat) variable
        end select
        if (iostat /= 0) then
          call error("failed to parse '" // setting // '=' // trim(string) // "' as a vector.")
        end if
      end if
    end subroutine

    subroutine string_print
      ! Print a config option
      write(string, '(a)') trim(setting)

      select type(variable)
        type is (real(wp))
          write(*, '(a,"= ",f10.5,",",/,22x,f10.5,",",/,22x,f10.5)') string(1:20), variable
      end select
    end subroutine
  end subroutine

  impure subroutine config_next(unit, iostat, string)
    ! This subroutine searches an initialization file for the next section and returns it.
    integer,            intent(in)  :: unit
    integer,            intent(out) :: iostat
    character(len=132)              :: string

    ! File and variable initialization
    string = ""
    iostat = 0

    ! Loop until a section is found or an error occurs
    do while (iostat == 0)
      ! Exit if the correct section is found
      if (string(1:1) == '[') then
        exit
      end if

      ! Read another line from the file
      call config_readline(unit, iostat, string)
    end do
  end subroutine

!subroutine structure_allocate(unit, conductors, superconductors, ferromagnets)
!  integer,              intent(in)  :: unit
!  integer, allocatable, intent(out) :: conductors(:)   
!  integer, allocatable, intent(out) :: superconductors(:)
!  integer, allocatable, intent(out) :: ferromagnets(:)
!
!  character(len=132)                :: string = ""
!  integer                           :: iostat = 0
!  integer                           :: c=0, s=0, f=0
!
!  ! Count the number of materials
!  do while (iostat == 0)
!    select case (adjustl(string))
!      case ('[conductor]')
!        c = c + 1
!      case ('[superconductor]')
!        s = s + 1
!      case ('[ferromagnet]')
!        f = f + 1
!    end select
!
!    read(unit, '(A)', iostat=iostat) string
!  end do
!  rewind(unit)
!
!  ! Allocate space for the materials
!  allocate(conductors(c), superconductors(s), ferromagnets(f))
!end subroutine

!subroutine structure_connect(unit, structure, conductors, superconductors, ferromagnets)
!  integer, pointer,     intent(out) :: structure
!
!  integer, pointer :: this, next
!
!  this%material_b => next
!  next%material_a => this
!end subroutine

!recursive subroutine update(structure)
!  class(material) :: structure
!
!  call structure%update
!  if (associated(structure%matrial_b))
!    call update(structure%material_b)
!    call structure%update
!  end if
!end subroutine
end module

