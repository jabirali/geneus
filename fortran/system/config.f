! This module define the interface 'config', which can be used to read and parse configuration settings from file or command line.
! TODO: This file is still incomplete. It should be finished, and then replace the option(...) routines used previously.
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
            call statement_substitute(statement % value)
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
    end subroutine

    subroutine statement_substitute(statement)
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
    end subroutine
  end function

  pure function config_count(name) result(count)
    character(len=132), intent(in) :: name
    integer                        :: count

    count = 0
  end function

  pure function config_find(name, number) result(count)
    character(len=132), intent(in) :: name
    integer,            intent(in) :: number
    integer                        :: count

    count = 0
  end function
end module

