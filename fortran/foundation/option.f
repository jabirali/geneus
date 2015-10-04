! This modules define the interface 'option', which can be used to read and parse command line arguments.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-10
! Updated: 2015-10-04

module mod_option
  use mod_math, only: wp
  implicit none
  private

  ! Public interface
  public option, print_option

  ! Generic procedures
  interface option
    module procedure option_logical, option_integer, option_real, option_reals, option_string
  end interface
contains
  subroutine option_integer(variable, option)
    ! Reads a command line option on the form option=value, where value is an integer.
    ! Note that 'variable' is only updated if the option is found, meaning that it should
    ! should be initialized to a sensible default value before this subroutine is called.
    integer,            intent(inout) :: variable
    character(len= * ), intent(in   ) :: option
    character(len=128)                :: string
    character(len= 22)                :: output
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
    real(kind=wp),      intent(inout) :: variable
    character(len= * ), intent(in   ) :: option
    character(len=128)                :: string
    character(len= 22)                :: output
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

  subroutine option_reals(array, option)
    ! Reads a command line option on the form option=value, where value is a real vector.
    ! Note that 'variable' is only updated if the option is found, meaning that it should
    ! should be initialized to a sensible default value before this subroutine is called.
    real(kind=wp),      intent(inout) :: array(3)
    character(len= * ), intent(in   ) :: option
    character(len=128)                :: string
    character(len= 22)                :: output
    integer                           :: n
    
    do n = 1,command_argument_count()
      ! Iterate over all command line arguments
      call get_command_argument(n,string)

      ! If this is the argument we were looking for, update the output variable
      if ( string(1:len(option)+1)  == option // '=' ) then
        read( string(len(option)+2:len(string)), * ) array(:)
      end if
    end do

    ! Write the results to standard out for verification purposes
    output = option
    write(*,'(a,a,f10.5,/,26x,f10.5,/,26x,f10.5)') ' :: ', output, array(:)
  end subroutine

  subroutine option_logical(variable, option)
    ! Reads a command line option on the form option=value, where value is a boolean.
    ! Note that 'variable' is only updated if the option is found, meaning that it should
    ! should be initialized to a sensible default value before this subroutine is called.
    logical,            intent(inout) :: variable
    character(len= * ), intent(in   ) :: option
    character(len=128)                :: string
    character(len= 22)                :: output
    integer                           :: n
    
    do n = 1,command_argument_count()
      ! Iterate over all command line arguments
      call get_command_argument(n,string)

      ! If this is the argument we were looking for, update the output variable
      if ( string(1:len(option)+1)  == option // '=' ) then
        read( string(len(option)+2:len(string)), '(l1)' ) variable
      end if
    end do

    ! Write the results to standard out for verification purposes
    output = option
    write(*,'(a,a,l10)') ' :: ', output, variable
  end subroutine

  subroutine option_string(variable, option)
    ! Reads a command line option on the form option=value, where value is a string.
    ! Note that 'variable' is only updated if the option is found, meaning that it should
    ! should be initialized to a sensible default value before this subroutine is called.
    character(len= * ), intent(inout) :: variable
    character(len= * ), intent(in   ) :: option
    character(len=128)                :: string
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
    write(*,'(a,a,1x,a,a,a)') ' :: ', option, '"', trim(variable), '"'
  end subroutine

  subroutine print_option
    ! Prints out a header for runtime parameters.
    write(*,'(a)') '╒═══════════════════════════════════╕'
    write(*,'(a)') '│        RUNTIME  PARAMETERS        │'
    write(*,'(a)') '╘═══════════════════════════════════╛'
  end subroutine
end module 
