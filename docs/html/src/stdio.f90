!> Author:   Jabir Ali Ouassou
!> Category: System
!>
!> This file renames the ISO input/output units to the standard UNIX names,
!> defines the ANSI escape codes for colored output, and defines a number of
!> auxiliary subroutines for e.g. writing out error and warning messages.

module stdio_m
    use, intrinsic :: iso_fortran_env
    public

    ! Declare standard input/output units
    integer :: stdin  = input_unit
    integer :: stdout = output_unit
    integer :: stderr = error_unit

    ! Define escape codes for terminal colors
    character(*), parameter :: color_none   = '[00m'
    character(*), parameter :: color_red    = '[31m'
    character(*), parameter :: color_green  = '[32m'
    character(*), parameter :: color_yellow = '[33m'
    character(*), parameter :: color_blue   = '[34m'
    character(*), parameter :: color_purple = '[35m'
    character(*), parameter :: color_cyan   = '[36m'
    character(*), parameter :: color_white  = '[37m'

    ! Declare public interfaces
    interface dump
        module procedure dump_arrays, dump_scalar
    end interface
contains
    subroutine message(msg)
    !!  Provides a way to report a status message.
        character(*), intent(in) :: msg

        write (stderr, '(a)') color_green // ' >> INFO: ' // color_none // msg
        flush (stderr)
    end subroutine

    subroutine warning(msg)
    !!  Provides a way to report a warning message.
        character(*), intent(in) :: msg

        write (stderr, '(a)') color_yellow // ' >> WARNING: ' // color_none // msg
        flush (stderr)
    end subroutine

    subroutine error(msg)
    !!  Provides a way to report an error message and halt the program.
        character(*), intent(in) :: msg

        write (stderr, '(a)') color_red // ' >> ERROR: ' // color_none // msg
        flush (stderr)
        stop
    end subroutine

    subroutine status_head(title)
    !!  Provides a way to write boxed status messages to standard out; in
    !!  particular, this routine writes out a boxed title with a timestamp.
        character(len=*), intent(in) :: title

        character(len=33) :: title_
        real              :: time
        integer           :: hh, mm, ss

        ! Calculate the current time
        call cpu_time(time)
        hh = int(time/3600.0)
        mm = int(mod(time, 3600.0)/60.0)
        ss = int(mod(time, 60.0))

        ! adjust the provided title
        title_ = ''
        title_((len(title_) - len(title) + 1)/2:) = title

        ! Write out the boxed header
        flush (unit=stdout)
        write (stdout, *)
        write (stdout, '(a)') &
            '╒═══════════════════════════════════╕'
        write (stdout, '(a)') &
            '│ ' // title_ // ' │'
        write (stdout, '(a)') &
            '├───────────────────────────────────┤'
        write (stdout, '(a,3x,a,7x,i3.2,a,i2.2,a,i2.2,3x,a)') &
            '│', 'Elapsed time:', hh, ':', mm, ':', ss, '│'
    end subroutine

    subroutine status_body(title, value)
    !!  Provides a way to write boxed status messages to standard out; in
    !!  particular, this routine writes out the name and value of a variable.
        character(len=*), intent(in) :: title
        class(*),         intent(in) :: value

        character(len=20) :: title_

        ! Adjust the provided title
        title_ = trim(title) // ':'

        ! Print out the title and value
        ! TODO: Use a kind selector when that is supported.
        select type (value)
            type is (integer)
                ! Standard integers
                write (stdout, '(a,3x,a,i10  ,2x,a)') '│', title_, value, '│'
            type is (real)
                ! Single-precision reals
                write (stdout, '(a,3x,a,f10.8,2x,a)') '│', title_, value, '│'
            type is (double precision)
                ! Double-precision reals
                write (stdout, '(a,3x,a,g10.8,2x,a)') '│', title_, value, '│'
        end select
    end subroutine

    subroutine status_foot()
    !!  Provides a way to write boxed status messages to standard out; in
    !!  particular, this routine writes out the bottom edge of such a box.

        ! Write out the boxed footer
        write (stdout, '(a)') &
            '╘═══════════════════════════════════╛'

        ! Flush the information to standard out
        flush (unit=stdout)
    end subroutine

    subroutine status_box(title)
    !!  Provides a way to write boxed status messages to standard out.
        character(len=*), intent(in) :: title

        character(len=33) :: title_

        ! Adjust the provided title
        title_ = ''
        title_((len(title_) - len(title) + 1)/2:) = title

        ! Write out the boxed message
        write (stdout, *)
        write (stdout, '(a)') &
            '╒═══════════════════════════════════╕'
        write (stdout, '(a)') &
            '│ ' // title_ // ' │'
        write (stdout, '(a)') &
            '╘═══════════════════════════════════╛'
    end subroutine

    function input(file) result(unit)
    !!  Open an input file for reading.
        character(len=*), intent(in) :: file
        integer                      :: unit

        integer :: iostat

        ! Open the output file
        open (newunit=unit, file=file, iostat=iostat, action='read', status='old')
        if (iostat /= 0) then
            call error('Failed to open input file "' // file // '"!')
        end if
    end function

    function output(file) result(unit)
    !!  Open an output file for writing.
        character(len=*), intent(in) :: file
        integer                      :: unit

        integer :: iostat

        ! Open the output file
        open (newunit=unit, file=file, iostat=iostat, action='write', status='replace')
        if (iostat /= 0) then
            call error('Failed to open output file "' // file // '"!')
        end if
    end function

    subroutine dump_arrays(filename, arrays, header)
    !!  Dump numerical arrays to an output file.
        use :: iso_fortran_env

        character(len=*),               intent(in) :: filename
        real(real64),     dimension(:), intent(in) :: arrays
        character(len=*), dimension(:), intent(in) :: header

        real(real64), dimension(size(arrays)/size(header), size(header)) :: matrix
        integer :: unit
        integer :: n

        ! Reshape the data
        matrix = reshape(arrays, shape(matrix))

        ! Open the output file
        unit = output(filename)

        ! Write the header line
        write (unit, '(*(a20,:,"        "))') '# ' // header(1), header(2:)

        ! Loop over the matrix rows
        do n = 1, size(matrix, 1)
            ! Write the matrix column to file
            write (unit, '(*(es20.12e3,:,"        "))') matrix(n, :)
        end do

        ! Close the output file
        close (unit=unit)
    end subroutine

    subroutine dump_scalar(filename, scalar)
    !!  Dump a numerical result to an output file.
        use :: iso_fortran_env

        character(len=*), intent(in) :: filename
        real(real64),     intent(in) :: scalar

        character(len=2048) :: str
        integer :: unit
        integer :: n

        ! Open the output file
        unit = output(filename)

        ! Write out command-line arguments
        do n = 1, command_argument_count()
            call get_command_argument(n, str)
            write (unit, '(a,"        ")', advance='no') trim(str)
        end do

        ! Write the scalar value to file
        write (unit, '(*(es20.12e3,:,"        "))') scalar

        ! Close the output file
        close (unit=unit)
    end subroutine
end module
