! This file defines a module containing the machine size of single-precision, double-precision, and quadruple-precision
! floating point numbers; to declare the floating point precision of a variable, use real(sp), real(dp), or real(qp) as
! the type of the variable. It also renames the ISO input/output units to the standard UNIX names.  Finally, the module
! defines a set of subroutines with the common interface 'option' to read and parse command line arguments to a program.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-10
! Updated: 2015-07-25

module mod_system
  use, intrinsic :: iso_fortran_env

  ! Declare standard input/output units
  integer,      parameter :: stdin   = input_unit
  integer,      parameter :: stdout  = output_unit
  integer,      parameter :: stderr  = error_unit

  ! Declare floating-point precisions
  integer,      parameter :: sp      = REAL32
  integer,      parameter :: dp      = REAL64
  integer,      parameter :: qp      = REAL128

  ! Define comm on mathematical constants
  real(dp),     parameter :: inf     = huge(1.0_dp)
  real(dp),     parameter :: pi      = atan(1.0_dp)*4.0_dp
  complex(dp),  parameter :: i       = (0.0_dp,1.0_dp)

  ! Define escape codes for terminal colors
  character(*), parameter :: color_none   = '[0m'
  character(*), parameter :: color_red    = '[31;1m'
  character(*), parameter :: color_green  = '[32;1m'
  character(*), parameter :: color_yellow = '[33;1m'
  character(*), parameter :: color_blue   = '[34;1m'
  character(*), parameter :: color_purple = '[35;1m'
  character(*), parameter :: color_cyan   = '[36;1m'
  character(*), parameter :: color_white  = '[37;1m'

  ! Define an interface for obtaining command line arguments
  interface option
    module procedure option_header, option_integer, option_real, option_string
  end interface
contains
  subroutine option_header
    ! If the subroutine 'option' is run without arguments, this prints out a header.
    write(*,'(a)') '╒═══════════════════════════════════╕'
    write(*,'(a)') '│        RUNTIME  PARAMETERS        │'
    write(*,'(a)') '╘═══════════════════════════════════╛'
  end subroutine

  subroutine option_integer(variable, option)
    ! Reads a command line option on the form option=value, where value is an integer.
    ! Note that 'variable' is only updated if the option is found, meaning that it should
    ! should be initialized to a sensible default value before this subroutine is called.
    integer,            intent(inout) :: variable
    character(len= * ), intent(in   ) :: option
    character(len=128)                :: string
    character(len= 20)                :: output
    integer                           :: n
    
    do n = 1,command_argument_count()
      ! Iterate over all command line arguments
      call get_command_argument(n,string)

      ! If this is the argument we were looking for, update the output variable
      if ( string(1:len(option)+1)  == option // '=' ) then
        read( string(len(option)+2:len(string)), '(i10)' ) variable
      end if
    end do

    ! Write the results to standard out for verification purposes
    output = option
    write(*,'(a,a,i10)') ' :: ', output, variable
  end subroutine

  subroutine option_real(variable, option)
    ! Reads a command line option on the form option=value, where value is a real number.
    ! Note that 'variable' is only updated if the option is found, meaning that it should
    ! should be initialized to a sensible default value before this subroutine is called.
    real(dp),           intent(inout) :: variable
    character(len= * ), intent(in   ) :: option
    character(len=128)                :: string
    character(len= 20)                :: output
    integer                           :: n
    
    do n = 1,command_argument_count()
      ! Iterate over all command line arguments
      call get_command_argument(n,string)

      ! If this is the argument we were looking for, update the output variable
      if ( string(1:len(option)+1)  == option // '=' ) then
        read( string(len(option)+2:len(string)), '(g24.0)' ) variable
      end if
    end do

    ! Write the results to standard out for verification purposes
    output = option
    write(*,'(a,a,f10.5)') ' :: ', output, variable
  end subroutine

  subroutine option_string(variable, option)
    ! Reads a command line option on the form option=value, where value is a string.
    ! Note that 'variable' is only updated if the option is found, meaning that it should
    ! should be initialized to a sensible default value before this subroutine is called.
    character(len= * ), intent(inout) :: variable
    character(len= * ), intent(in   ) :: option
    character(len=128)                :: string
    character(len= 20)                :: output
    integer                           :: n
    
    do n = 1,command_argument_count()
      ! Iterate over all command line arguments
      call get_command_argument(n,string)

      ! If this is the argument we were looking for, update the output variable
      if ( string(1:len(option)+1)  == option // '=' ) then
        read( string(len(option)+2:len(string)), '(a)' ) variable
      end if
    end do

    ! Write the results to standard out for verification purposes
    output = option
    write(*,'(a,a,a)') ' :: ', output, variable
  end subroutine
end module 
