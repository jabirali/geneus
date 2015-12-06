! This file renames the ISO input/output units to the standard UNIX names, and defines the ANSI escape codes for colored output.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-10
! Updated: 2015-10-04

module mod_stdio
  use iso_fortran_env
  implicit none
  public

  ! Declare standard input/output units
  integer,      parameter :: stdin   = input_unit
  integer,      parameter :: stdout  = output_unit
  integer,      parameter :: stderr  = error_unit

  ! Define escape codes for terminal colors
  character(*), parameter :: color_none   = '[0m'
  character(*), parameter :: color_red    = '[31m'
  character(*), parameter :: color_green  = '[32m'
  character(*), parameter :: color_yellow = '[33m'
  character(*), parameter :: color_blue   = '[34m'
  character(*), parameter :: color_purple = '[35m'
  character(*), parameter :: color_cyan   = '[36m'
  character(*), parameter :: color_white  = '[37m'
contains
  impure subroutine message(msg)
    ! This subroutine provides a way to report a status message.
    character(*), intent(in) :: msg

    write(stderr,'(a)') color_green  // '>> INFO: '    // color_none // msg
  end subroutine

  impure subroutine warning(msg)
    ! This subroutine provides a way to report a warning message.
    character(*), intent(in) :: msg

    write(stderr,'(a)') color_yellow // '>> WARNING: ' // color_none // msg
  end subroutine

  impure subroutine error(msg)
    ! This subroutine provides a way to report an error message and halt the program.
    character(*), intent(in) :: msg

    write(stderr,'(a)') color_red    // '>> ERROR: '   // color_none // msg
    stop
  end subroutine
end module 
