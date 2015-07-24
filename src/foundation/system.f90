! This file defines a module containing the machine size of single-precision, double-precision, and quadruple-precision
! floating point numbers; to declare the floating point precision of a variable, use real(sp), real(dp), or real(qp) as
! the type of the variable.  Furthermore, the module defines a set of subroutines with the common interface 'option' to
! read and parse command line arguments to a program. Finally, the module defines the two subroutines 'print_error' and 
! and 'print_info' for writing error messages and other information to standard out.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-10
! Updated: 2015-07-24

module mod_system
  use, intrinsic :: iso_fortran_env

  ! Declare floating-point precisions
  integer,  parameter :: sp  = REAL32
  integer,  parameter :: dp  = REAL64
  integer,  parameter :: qp  = REAL128

  ! Define common mathematical constants
  real(dp), parameter :: inf = huge(1.0_dp)
  real(dp), parameter :: pi  = atan(1.0_dp)*4.0_dp

  ! Define an interface for obtaining command line arguments
  interface option
    module procedure option_header, option_integer, option_real, option_string
  end interface
contains
  subroutine option_header
    ! If the subroutine 'option' is run without arguments, this prints an option header.
    write(*,'(a)') '╒═══════════════════════════════════╕'
    write(*,'(a)') '│        RUNTIME  PARAMETERS        │'
    write(*,'(a)') '╘═══════════════════════════════════╛'
  end subroutine

  subroutine option_integer(variable, option)
    ! Reads a command line option on the form --option=value, where value is an integer.
    ! Note that 'variable' is only updated if the option is found, meaning that it should
    ! should be initialized to a sensible default value before this subroutine is called.
    integer,            intent(inout) :: variable
    character(len= * ), intent(in   ) :: option
    character(len=128)                :: string
    character(len= 20)                :: output
    
    do n = 1,command_argument_count()
      ! Iterate over all command line arguments
      call get_command_argument(n,string)

      ! If this is the argument we were looking for, update the output variable
      if ( string(1:len(option)+3)  == '--' // option // '=' ) then
        read( string(len(option)+4:len(string)), '(i10)' ) variable
      end if
    end do

    ! Write the results to standard out for verification purposes
    output = option
    write(*,'(a,a,i10)') ' :: ', output, variable
  end subroutine

  subroutine option_real(variable, option)
    ! Reads a command line option on the form --option=value, where value is a real number.
    ! Note that 'variable' is only updated if the option is found, meaning that it should
    ! should be initialized to a sensible default value before this subroutine is called.
    real(dp),           intent(inout) :: variable
    character(len= * ), intent(in   ) :: option
    character(len=128)                :: string
    character(len= 20)                :: output
    
    do n = 1,command_argument_count()
      ! Iterate over all command line arguments
      call get_command_argument(n,string)

      ! If this is the argument we were looking for, update the output variable
      if ( string(1:len(option)+3)  == '--' // option // '=' ) then
        read( string(len(option)+4:len(string)), '(g24.0)' ) variable
      end if
    end do

    ! Write the results to standard out for verification purposes
    output = option
    write(*,'(a,a,f10.5)') ' :: ', output, variable
  end subroutine

  subroutine option_string(variable, option)
    ! Reads a command line option on the form --option=value, where value is a string.
    ! Note that 'variable' is only updated if the option is found, meaning that it should
    ! should be initialized to a sensible default value before this subroutine is called.
    character(len= * ), intent(inout) :: variable
    character(len= * ), intent(in   ) :: option
    character(len=128)                :: string
    character(len= 20)                :: output
    
    do n = 1,command_argument_count()
      ! Iterate over all command line arguments
      call get_command_argument(n,string)

      ! If this is the argument we were looking for, update the output variable
      if ( string(1:len(option)+3)  == '--' // option // '=' ) then
        read( string(len(option)+4:len(string)), '(a)' ) variable
      end if
    end do

    ! Write the results to standard out for verification purposes
    output = option
    write(*,'(a,a,a)') ' :: ', output, variable
  end subroutine

  subroutine print_error(str)
    ! This subroutine can be used to write error messages to standard out.
    character(len=*), intent(in) :: str
    
    write(*,'(a1,a,a1,a,a)') achar(27), '[1;31mERROR:', achar(27), '[0m ', str
  end subroutine

  subroutine print_info(str)
    ! This subroutine can be used to write info messages to standard out.
    character(len=*), intent(in) :: str
    
    write(*,'(a1,a,a1,a,a)') achar(27), '[1;33mINFO: ', achar(27), '[0m ', str
  end subroutine
end module 
