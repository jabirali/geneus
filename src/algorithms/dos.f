!module mod_dos
!  use mod_system
!  use mod_material
!  use mod_multilayer
!  implicit none
!contains
!
!subroutine calculate_dos(material, iterations, unit)
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  class(material),  target, intent(inout) :: material
!  integer,         optional, intent(in   ) :: iterations
!  integer,         optional, intent(in   ) :: unit
!
!  class(material), pointer                :: p           => null()
!  integer                                  :: iterations_ =  0
!  integer                                  :: unit_       =  stdout
!  integer                                  :: n 
!
!  ! Handle optional arguments
!  if (present(iterations)) then
!    iterations_ = iterations
!  end if
!
!  if (present(unit)) then
!    unit_ = unit
!  end if
!
!  ! Update the state of the structure
!  if (iterations_ <= 0) then
!    ! If iterations is zero, update the provided material only
!    call material % update
!
!    ! Update the output file
!    ! call print_results
!  else
!    ! Update the entire multilayer structure if 'iterations' is positive
!    do n = 1,iterations
!      call print_progress(n)
!      call update_all(material)
!
!    ! Update the output file
!    ! call print_results
!    end do
!  end if
!
!  ! Iterate to the left end of the multilayer structure
!  p => material
!  do while (associated(p % material_a))
!    p => p % material_a
!  end do
!
!  ! Iterate to the right end of the structure, and write
!  ! out the density of states of each material on the way
!  n = 0
!  do while (associated(p))
!    call p % write_dos(unit_, dble(n), dble(n+1))
!    p => p % material_b
!    n =  n+1
!  end do
!contains
!  subroutine print_progress(n)
!    ! !!!!
!    real(sp) :: time
!    integer  :: n
!    
!    ! Determine how much CPU time has elapsed
!    call cpu_time(time)
!
!    ! Print the progress information to standard out
!    write(*,'(a)') '                                     '
!    write(*,'(a)') '╒═══════════════════════════════════╕'
!    write(*,'(a)') '│       PROGRESS  INFORMATION       │'
!    write(*,'(a)') '├───────────────────────────────────┤'
!    write(*,'(a,6x,a,i2.2,a,i2.2,7x,a)')                &
!      '│','Iteration:     ', n, ' / ', iterations_,    '│'
!    write(*,'(a,6x,a,i2.2,a,i2.2,a,i2.2,6x,a)')         &
!      '│','Elapsed time:  ',                            &
!      int(time/3600.0_sp),':',                          &
!      int(mod(time,3600.0_sp)/60.0_sp),':',             &
!      int(mod(time,60.0_sp)),                          '│'
!    write(*,'(a)') '╘═══════════════════════════════════╛'
!    
!  end subroutine 
!end subroutine
!
!
!end module
!
