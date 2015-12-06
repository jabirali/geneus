! This module define the interface 'config', which can be used to read and parse configuration settings from file or command line.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-12-03
! Updated: 2015-12-05

module mod_config 
  use mod_stdio, only: message, warning, error
  use mod_math,  only: wp
  implicit none
  !private

  type :: statement_t
    character(len=132) :: field
    character(len=132) :: value
  end type

  type :: section_t
    character(len=132) :: name
    integer            :: size
    type(statement_t)  :: statement(128)
  end type

  type :: config_t
    character(len=132)          :: file
    integer                     :: unit
    integer                     :: iostat
    integer                     :: size
    type(section_t)             :: section(128)
  end type

  interface config_t
    module procedure config_construct
  end interface

  interface config
    module procedure config_scalar, config_vector, config_next
  end interface

  !public config
contains



  ! NEW OBJECT ORIENTED APPROACH



  impure function config_construct(file) result(config)
    character(len= * ), intent(in) :: file
    type(config_t)                 :: config
    character(len=132)             :: string
    integer                        :: n, m, i, j, k

    ! Status message
    call message('loading config file "' // file // '"...')

    ! Initialize the data type
    config % file    = file
    config % unit    = 0
    config % iostat  = 0
    config % size    = 0

    ! Open the config file
    open(newunit = config % unit, file = config % file, iostat = config % iostat, action = 'read', status = 'old')
    if (config % iostat /= 0) then
      call error('failed to open configuration file!')
    end if

    ! Read the config file
    do while (config % iostat == 0)
      ! Read one line
      string = ''
      read(config % unit, '(a)', iostat = config % iostat) string

      ! Process line
      if (scan(string, '[]') /= 0) then
        ! If the line has the format [...], process it as a section
        call section_register
      else if (scan(string, '=') /= 0) then
        ! If the line has the format a = b, process it as a statement
        call statement_register
      else if (config % size > 0) then
        if (config % section(config%size) % size > 0) then
          ! If possible, append it to the previous statement
          call statement_append
        end if
      end if
    end do

    ! Close the config file
    close(unit = config % unit)

    ! Substitute command line arguments and check for duplicates
    do n = 1,config%size
      associate(section => config % section(n))
        do m = 1,section%size
          associate(statement => section % statement(m))
            call substitute(statement % value)
            do k = 1,m-1
              associate(current => section % statement(k))
                if (current % field == statement % field) then
                  call error('multiple definitions of field "' // trim(statement % field) // '" for ['// trim(section%name) //'].')
                end if
              end associate
            end do
          end associate
        end do
      end associate
    end do

    ! Status message
    call message('successfully parsed config file!')
  contains
    subroutine section_register
      ! Find the section delimiters
      i = scan(string, '[')
      j = scan(string, ']')
      if (i <= 0 .or. j <= i) then
        call error('failed to parse section "' // trim(adjustl(string)) // '".')
      end if

      ! Register a new section
      config  % size = config % size + 1
      associate(section => config % section(config % size))
        section % name = adjustl(string(i+1:j-1))
        section % size = 0
      end associate

      write(*,*) config % section(config % size) % name, config % size
    end subroutine

    subroutine statement_register
      ! Find the assignment delimiter
      i = scan(string, '=')
      if (string(1:i-1) == '' .or. string(i+1:) == '') then
        call error('failed to parse assignment "' // trim(adjustl(string)) // '".')
      end if

      ! Find the comment delimiter
      j = scan(string, '#')
      if (j == 0) then
        j = len(string)
      else if (j < i) then
        return
      end if

      ! Register a new statement
      associate(section   => config % section(config%size))
      section   % size  = section % size + 1
        associate(statement => section % statement(section%size))
          statement % field = adjustl(string(  1:i-1))
          statement % value = adjustl(string(i+1:j-1))
        end associate
      end associate

      write(*,*) config % section(config%size) % size
      write(*,*) config % section(config%size) % statement(config%section(config%size)%size) % field 
      write(*,*) config % section(config%size) % statement(config%section(config%size)%size) % value
    end subroutine

    subroutine statement_append
      ! Find the comment delimiter
      j = scan(string, '#')
      if (j == 0) then
        j = len(string)
      end if

      ! Append data to current statement
      associate(section   => config  % section(config%size))
        associate(statement => section % statement(section%size))
          statement % value = trim(statement % value) // trim(adjustl(string(1:j)))
        end associate
      end associate

      write(*,*) config % section(config%size) % statement(config%section(config%size)%size) % value
    end subroutine

    subroutine substitute(statement)
      character(len=132), intent(inout) :: statement
      integer                           :: length

      ! Find the command line argument marker
      i = scan(statement, '$')
      if (i == 0) then
        return
      end if

      ! Read the command line argument number
      j = 0
      read(statement(i+1:),*,iostat=config%iostat) j
      if (j <= 0) then
        call error('invalid command line argument "' // trim(statement(i:)) // '"!')
      end if

      ! Read the command line argument itself
      string = ''
      call get_command_argument(number=j, value=string, length=length)
      if (length <= 0) then
        call error('failed to obtain command line arguments!')
      end if

      ! Replace the statement with the argument
      statement = string

      write(*,*) statement
    end subroutine
  end function



  ! OLDER PROCEDURAL APPROACH




  impure subroutine config_readline(unit, iostat, output, reverse)
    ! This subroutine reads one non-empty line from an initialization file, strips it of comments and whitespace, and returns it.
    ! If an error occurs during reading, such as an end-of-file or read error, this will be reflected by a non-zero iostat value.
    integer,            intent(in)  :: unit
    integer,            intent(out) :: iostat
    logical, optional,  intent(in)  :: reverse
    character(len=132), intent(out) :: output
    character(len=132)              :: input
    character(len=132)              :: string
    integer                         :: conts
    integer                         :: n, m, i, j

    ! Variable initialization
    input  = ""
    output = ""
    iostat = 0
    conts  = 0
    m      = 1

    ! Process one non-empty line
    do while (output == "")
      ! Move backwards if 'reverse' is enabled
      if (present(reverse)) then
        if (reverse) then
          backspace(unit)
          backspace(unit)
        end if
      end if

      ! Read one line from the input file
      read(unit,'(a)',iostat=iostat) input
      if (iostat /= 0) then
        exit
      end if

      ! Skip lines without sections or statements
      if (scan(input,'[=]') == 0) then
        cycle
      end if

      ! Parse the current line
      n = 0
      do
        ! Make sure we stop at the end of the line
        n = n + 1
        if (n > len(input)) then
          exit
        end if

        ! Process the current character of the line
        select case(input(n:n))
          ! Skip space and tab characters
          case (' ', char(11))
            cycle

          ! '#' is interpreted as a start of a comment
          case ('#')
            exit

          ! '$n' is interpreted as a command line argument
          case ('$')
            ! Interpret the next character as the argument number
            n = n + 1
            read(input(n:n),*,iostat=iostat) i
            if (iostat /= 0) then
              call error("failed to parse config file: bad argument '$" // input(n:n) // "'.")
            end if

            ! Read the corresponding command line argument
            call get_command_argument(number=i, value=string, length=j)
            if (j <= 0) then
              call error("failed to parse config file: missing argument '$" // input(n:n) // "'.")
            end if

            ! Update the output string 
            output(m:m+j) = string
            m = m + j
            cycle

          ! ',' may require a line continuation
          case (',')
            ! Read another line if neccessary
            if (scan(input(n+1:),'0123456789') == 0) then
              n = 1
              conts = conts + 1
              read(unit,'(a)',iostat=iostat) input
              if (iostat /= 0) then
                call error("failed to parse config file: file ended with a comma.")
              end if
            end if

            ! Update the output string
            output(m:m) = ','
            m = m + 1
            cycle

          ! Other characters require no special handling
          case default
            output(m:m) = input(n:n)
            m = m + 1
            cycle
        end select
      end do
    end do

    ! Revert continuations by moving backwards
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
          call error("failed to parse config file: '" // setting // '=' // trim(string) // "' is not a valid scalar assignment.")
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
          call error("failed to parse config file: '" // setting // '=' // trim(string) // "' is not a valid vector assignment.")
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

