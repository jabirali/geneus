! This program calculates the density of states in an S/N/F/S Josephson Junction, assuming that both of
! the superconductors can be described as BCS superconductors with a set phase difference.  The program
! works by initially setting the phase difference to zero,   then bootstrapping the system by gradually
! increasing the phase difference to the provided value, and then continuing to update the state of the
! system until the Green's functions of the middle layers converge. The density of states as a function
! of position and energy is then written to an output file. Note that the output file may be visualized
! using the accompanying Gnuplot script plot/dos.plt.
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-08-04
! Updated: 2015-08-05

program dos_snfs
  use mod_hybrid
  implicit none

  !--------------------------------------------------------------------------------!
  !                            DECLARATION OF VARIABLES                            !
  !--------------------------------------------------------------------------------!

  ! Declare the materials in the hybrid structure
  type(superconductor)    :: s1
  type(conductor)         :: n
  type(ferromagnet)       :: f
  type(superconductor)    :: s2

  ! Declare the variables used by the program
  integer                 :: output
  real(dp),   allocatable :: energy_array(:)
  real(dp)                :: interfaces(5)
  real(dp)                :: phase_
  integer                 :: iteration

  ! Declare the parameters that can be modified at runtime
  character(len=32)       :: filename             = 'dos_snfs.dat'
  character(len=32)       :: structure            = 'SNFS'
  integer                 :: bootstraps           =  15
  integer                 :: energies             = 150      
  real(dp)                :: scattering           = 0.01_dp
  real(dp)                :: conductance          = 0.30_dp
  real(dp)                :: phase                = 0.00_dp
  real(dp)                :: s_length             = 1.00_dp
  real(dp)                :: s_coupling           = 0.20_dp
  logical                 :: n_enable             = .false.
  real(dp)                :: n_length             = 0.50_dp
  real(dp)                :: n_spinorbit_a        = 0.00_dp
  real(dp)                :: n_spinorbit_b        = 0.00_dp
  logical                 :: f_enable             = .false.
  real(dp)                :: f_length             = 0.50_dp
  real(dp)                :: f_exchange_x         = 0.00_dp
  real(dp)                :: f_exchange_y         = 0.00_dp
  real(dp)                :: f_exchange_z         = 0.00_dp



  !--------------------------------------------------------------------------------!
  !                           INPUT/OUTPUT PREPARATIONS                            !
  !--------------------------------------------------------------------------------!

  ! Process command line options
  call option
  call option(filename,             'filename')
  call option(structure,            'structure')

  ! Check the selected hybrid structure
  select case (trim(structure))
    case ('SNFS', 'SFNS')
      write(*,'(/,a,/)') ' :: '//color_green//trim(structure)//color_none//' structures are supported.'
      f_enable = .true.
      n_enable = .true.
    case ('SFS')
      write(*,'(/,a,/)') ' :: '//color_green//trim(structure)//color_none//' structures are supported.'
      f_enable = .true.
    case ('SNS')
      write(*,'(/,a,/)') ' :: '//color_green//trim(structure)//color_none//' structures are supported.'
      n_enable = .true.
    case default
      write(*,'(/,a,/)') ' :: '//color_red  //trim(structure)//color_none//' structures are not supported!'
      stop
  end select

  ! Process command line options
  call option(bootstraps,           'bootstraps')
  call option(energies,             'energies')
  call option(scattering,           'scattering')
  call option(conductance,          'conductance')
  call option(phase,                'phase')
  call option(s_length,             's.length')
  call option(s_coupling,           's.coupling')
  select case (trim(structure))
    case ('SNFS')
      call option(n_length,             'n.length')
      call option(n_spinorbit_a,        'n.spinorbit.a')
      call option(n_spinorbit_b,        'n.spinorbit.b')
      call option(f_length,             'f.length')
      call option(f_exchange_x,         'f.exchange.x')
      call option(f_exchange_y,         'f.exchange.y')
      call option(f_exchange_z,         'f.exchange.z')
    case ('SFNS')
      call option(f_length,             'f.length')
      call option(f_exchange_x,         'f.exchange.x')
      call option(f_exchange_y,         'f.exchange.y')
      call option(f_exchange_z,         'f.exchange.z')
      call option(n_length,             'n.length')
      call option(n_spinorbit_a,        'n.spinorbit.a')
      call option(n_spinorbit_b,        'n.spinorbit.b')
    case ('SNS')
      call option(n_length,             'n.length')
      call option(n_spinorbit_a,        'n.spinorbit.a')
      call option(n_spinorbit_b,        'n.spinorbit.b')
    case ('SFS')
      call option(f_length,             'f.length')
      call option(f_exchange_x,         'f.exchange.x')
      call option(f_exchange_y,         'f.exchange.y')
      call option(f_exchange_z,         'f.exchange.z')
  end select

  ! Open the output file
  open(newunit=output, file=filename)




  !--------------------------------------------------------------------------------!
  !                           VARIABLE INITIALIZATION                              !
  !--------------------------------------------------------------------------------!

  ! Allocate and initialize the energy array
  allocate(energy_array(energies))
  call energy_range(energy_array)

  ! Calculate the locations of the interfaces
  interfaces(1) = 0.0_dp
  !do iteration=2,len(structure)
  !  select case (structure(iteration))
  !    case ('S')
  !      interfaces(iteration) = interfaces(iteration-1) + s_length
  !    case ('N')
  !      interfaces(iteration) = interfaces(iteration-1) + n_length
  !    case ('F')
  !      interfaces(iteration) = interfaces(iteration-1) + f_length
  !  end select
  !end do

  ! Initialize the material layers
  s1 = superconductor (energy_array, scattering = scattering, thouless = 1/s_length**2, coupling  = s_coupling)
  s2 = superconductor (energy_array, scattering = scattering, thouless = 1/s_length**2, coupling  = s_coupling)
  n  = conductor      (energy_array, scattering = scattering, thouless = 1/n_length**2, &
                       spinorbit = spinorbit_xy(alpha = n_spinorbit_a, beta = n_spinorbit_b))
  f  = ferromagnet    (energy_array, scattering = scattering, thouless = 1/f_length**2, &
                       exchange  = [f_exchange_x, f_exchange_y, f_exchange_z])

  ! Connect the material layers
  select case (trim(structure))
    case('SNFS')
      call connect(s1, n, conductance, conductance)
      call connect(n,  f, conductance, conductance)
      call connect(f, s2, conductance, conductance)
    case('SFNS')
      call connect(s1, f, conductance, conductance)
      call connect(f,  n, conductance, conductance)
      call connect(n, s2, conductance, conductance)
    case('SFS')
      call connect(s1, f, conductance, conductance)
      call connect(f, s2, conductance, conductance)
    case('SNS')
      call connect(s1, n, conductance, conductance)
      call connect(n, s2, conductance, conductance)
  end select



  !--------------------------------------------------------------------------------!
  !                           BOOTSTRAPPING PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Write the initial density of states to file
  call write_dos(output)

  ! Bootstrap the system at zero phase difference
  if (n_enable .and. f_enable) then
    do iteration=1,3
      ! Status information
      call print_init

      ! Update the state of the system
      call f%update
      call n%update

      ! Write the density of states to file
      call write_dos(output)
    end do
  end if

  ! Bootstrap the system at gradually increasing phase difference
  do iteration=1,bootstraps
    ! Calculate the current phase difference
    phase_ = (phase*iteration)/bootstraps

    ! Status information
    call print_boot

    ! Reinitialize the superconductors
    call s1%init( exp((-i*pi/2)*phase_) )
    call s2%init( exp((+i*pi/2)*phase_) )

    ! Update the state of the system
    if (f_enable) then
      call f%update
    end if

    if (n_enable) then
      call n%update
    end if

    ! Write the density of states to file
    call write_dos(output)
  end do


  !--------------------------------------------------------------------------------!
  !                            CONVERGENCE PROCEDURE                               !
  !--------------------------------------------------------------------------------!

  ! Iterate until convergence
  iteration = 0
  do while (n%difference > 10*n%tolerance .or. f%difference > 10*f%tolerance)
    iteration = iteration + 1

    ! Status information
    call print_main

    ! Update the state of the system
    if (f_enable) then
      call f%update
    end if

    if (n_enable) then
      call n%update
    end if

    ! Write the density of states to file
    call write_dos(output)
  end do



  !--------------------------------------------------------------------------------!
  !                              CLEANUP PROCEDURE                                 !
  !--------------------------------------------------------------------------------!

  ! Close the output file
  close(unit=output)

  ! Deallocate memory
  deallocate(energy_array)

contains

  !--------------------------------------------------------------------------------!
  !                           INPUT/OUTPUT PROCEDURES                              !
  !--------------------------------------------------------------------------------!

  subroutine print_init
    ! Determine how much CPU time has elapsed
    real(sp) :: time
    call cpu_time(time)

    ! Print the progress information to standard out
    write(*,'(a)') '                                     '
    write(*,'(a)') '╒═══════════════════════════════════╕'
    write(*,'(a)') '│          INITIALIZATION           │'
    write(*,'(a)') '├───────────────────────────────────┤'
    write(*,'(a,5x,a,i8,4x,a)')                         &
      '│','Iteration:        ',iteration,'│'
    write(*,'(a,5x,a,f8.6,4x,a)')                       &
      '│','Phase difference: ',phase_+1e-12,'│'
    write(*,'(a,5x,a,i2.2,a,i2.2,a,i2.2,4x,a)')         &
      '│','Elapsed time:     ',                         &
      int(time/3600.0_sp),':',                          &
      int(mod(time,3600.0_sp)/60.0_sp),':',             &
      int(mod(time,60.0_sp)),                          '│'
    write(*,'(a)') '╘═══════════════════════════════════╛'
  end subroutine

  subroutine print_boot
    ! Determine how much CPU time has elapsed
    real(sp) :: time
    call cpu_time(time)

    ! Print the progress information to standard out
    write(*,'(a)') '                                     '
    write(*,'(a)') '╒═══════════════════════════════════╕'
    write(*,'(a)') '│          BOOTSTRAPPING            │'
    write(*,'(a)') '├───────────────────────────────────┤'
    write(*,'(a,5x,a,i8,4x,a)')                         &
      '│','Iteration:        ',iteration,'│'
    write(*,'(a,5x,a,f8.6,4x,a)')                       &
      '│','Phase difference: ',phase_+1e-12,'│'
    write(*,'(a,5x,a,i2.2,a,i2.2,a,i2.2,4x,a)')         &
      '│','Elapsed time:     ',                         &
      int(time/3600.0_sp),':',                          &
      int(mod(time,3600.0_sp)/60.0_sp),':',             &
      int(mod(time,60.0_sp)),                          '│'
    write(*,'(a)') '╘═══════════════════════════════════╛'
  end subroutine

  subroutine print_main
    ! Determine how much CPU time has elapsed
    real(sp) :: time
    call cpu_time(time)

    ! Print the progress information to standard out
    write(*,'(a)') '                                     '
    write(*,'(a)') '╒═══════════════════════════════════╕'
    write(*,'(a)') '│           CONVERGENCE             │'
    write(*,'(a)') '├───────────────────────────────────┤'
    write(*,'(a,5x,a,i8,4x,a)')                         &
      '│','Iteration:        ',iteration,'│'
    write(*,'(a,5x,a,f8.6,4x,a)')                       &
      '│','Maximum change:   ',max(n%difference,f%difference),'│'
    write(*,'(a,5x,a,i2.2,a,i2.2,a,i2.2,4x,a)')         &
      '│','Elapsed time:     ',                         &
      int(time/3600.0_sp),':',                          &
      int(mod(time,3600.0_sp)/60.0_sp),':',             &
      int(mod(time,60.0_sp)),                          '│'
    write(*,'(a)') '╘═══════════════════════════════════╛'
  end subroutine

  subroutine write_dos(output)
    integer, intent(in) :: output

    rewind(unit=output)
    select case (trim(structure))
      case ('SNFS')
        call s1 % write_dos(output, interfaces(1), interfaces(2))
        call n  % write_dos(output, interfaces(2), interfaces(3))
        call f  % write_dos(output, interfaces(3), interfaces(4))
        call s2 % write_dos(output, interfaces(4), interfaces(5))
      case ('SFNS')
        call s1 % write_dos(output, interfaces(1), interfaces(2))
        call f  % write_dos(output, interfaces(3), interfaces(4))
        call n  % write_dos(output, interfaces(2), interfaces(3))
        call s2 % write_dos(output, interfaces(4), interfaces(5))
      case ('SFS')
        call s1 % write_dos(output, interfaces(1), interfaces(2))
        call f  % write_dos(output, interfaces(3), interfaces(4))
        call s2 % write_dos(output, interfaces(4), interfaces(5))
      case ('SNS')
        call s1 % write_dos(output, interfaces(1), interfaces(2))
        call n  % write_dos(output, interfaces(2), interfaces(3))
        call s2 % write_dos(output, interfaces(4), interfaces(5))
    end select
    flush(unit=output)
  end subroutine
end program
