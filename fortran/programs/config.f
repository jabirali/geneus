module mod_config 
  implicit none
  public config_file_section, config_file_setting
  private
contains
  impure subroutine config_file_readline(unit, iostat, output)
    ! This subroutine reads one non-empty line from an initialization file, strips it of comments and whitespace, and returns it.
    ! If an error occurs during reading, such as an end-of-file or read error, this will be reflected by a non-zero iostat value.
    integer,            intent(in)  :: unit
    integer,            intent(out) :: iostat
    character(len=132), intent(out) :: output
    character(len=132)              :: input
    integer                         :: n, m

    ! Variable initialization
    m      = 1
    iostat = 0
    output = ""
    input  = ""

    ! Process one non-empty line before returning
    do while (output == "")
      ! Read one line from file
      read(unit,'(a)',iostat=iostat) input
      if (iostat /= 0) then
        exit
      end if

      ! Process the input line
      do n=1,len(input)
        ! Ignore comments
        if (input(n:n) == '#') then
          exit
        end if

        ! Ignore whitespace
        if (input(n:n) == ' ') then
          cycle
        end if

        ! Copy characters to output string
        output(m:m) = input(n:n) 
        m = m + 1
      end do
    end do
  end subroutine

  impure subroutine config_file_section(unit, iostat, section)
    ! This subroutine searches an initialization file for a given section, and returns an iostat:
    ! iostat = 0 if section found, iostat < 0 if end of file reached, iostat > 0 for misc errors.
    integer,            intent(in)  :: unit
    integer,            intent(out) :: iostat
    character(len= * ), intent(in)  :: section
    character(len=132)              :: string

    ! File and variable initialization
    rewind(unit)
    string = ""
    iostat = 0

    ! Loop until the section is found or an error occurs
    do while (iostat == 0)
      ! Exit if the correct section is found
      if (string == '['//section//']') then
        exit
      end if

      ! Read another line from the file
      call config_file_readline(unit, iostat, string)
    end do
  end subroutine

  impure subroutine config_file_setting(unit, iostat, setting, value)
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
      call config_file_readline(unit, iostat, string)
    end do

    ! Rewind until the beginning of the section
    string = ""
    do while (string(1:1) /= '[')
      backspace(unit)
      backspace(unit)
      read(unit,'(1a)') string
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

!program config
!  use mod_config
!  implicit none
!
!  integer            :: unit
!  integer            :: iostat = 0
!  character(len=132) :: string = ""
!
!  integer :: c, s, f
!
!  integer, pointer     :: structure
!  integer, allocatable :: conductors(:)   
!  integer, allocatable :: superconductors(:)
!  integer, allocatable :: ferromagnets(:)
!
!  ! Open the configuration file
!  open(newunit=unit, file="config.ini")
!
!  ! Allocate memory
!  call structure_allocate(unit, conductors, superconductors, ferromagnets)
!
!  ! Close the configuration file
!  close(unit)
!
!  write(*,*) size(conductors), size(superconductors), size(ferromagnets)

program config
  use mod_config
  implicit none

  integer :: iostat = 0
  integer :: unit   = 0
  character(len=132) :: string = ""

  open(newunit=unit,file='config.ini')
  call config_file_section(unit, iostat, 'general')
  write(*,*) iostat
  call config_file_setting(unit, iostat, 'magnetization',   string)
  write(*,*) string
  write(*,*) iostat
  call config_file_setting(unit, iostat, 'ferromagnets',    string)
  write(*,*) string
  write(*,*) iostat
  call config_file_setting(unit, iostat, 'superconductors', string)
  write(*,*) string
  write(*,*) iostat
  call config_file_setting(unit, iostat, 'magnetization',   string)
  write(*,*) string
  write(*,*) iostat
  call config_file_setting(unit, iostat, 'superconductors', string)
  write(*,*) string
  write(*,*) iostat
end program

