! This program calculates the density of states in an S/N/F trilayer, where the superconducting layer may
! be treated either self-consistently or not.  The density of states as a function of position and energy
! is then written to an output file.   Note that the output file may be visualized using the accompanying
! Gnuplot script 'plot/dos.plt'. 
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-08-05
! Updated: 2015-08-05

program dos_snf
  use mod_hybrid
  implicit none

  !--------------------------------------------------------------------------------!
  !                            DECLARATION OF VARIABLES                            !
  !--------------------------------------------------------------------------------!

  ! Declare the materials in the hybrid structure
  type(superconductor)    :: s
  type(conductor)         :: n
  type(ferromagnet)       :: f

  ! Declare the variables used by the program
  integer                 :: output
  real(dp),   allocatable :: energy_array(:)
  real(dp)                :: interfaces(4)
  integer                 :: iteration

  ! Declare the parameters that can be modified at runtime
  character(len=64)       :: filename             = 'dos_snf.dat'
  logical                 :: selfconsistent       = .false.
  integer                 :: energies             
  real(dp)                :: scattering           = 0.01_dp
  real(dp)                :: conductance          = 0.30_dp
  real(dp)                :: s_length             = 1.00_dp
  real(dp)                :: s_coupling           = 0.20_dp
  real(dp)                :: n_length             = 0.50_dp
  real(dp)                :: n_spinorbit_a        = 0.00_dp
  real(dp)                :: n_spinorbit_b        = 0.00_dp
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
  call option(selfconsistent,       'selfconsistent')

  if (selfconsistent) then
    energies = 600
  else
    energies = 150
  end if

  call option(energies,             'energies')
  call option(scattering,           'scattering')
  call option(conductance,          'conductance')
  call option(s_length,             's.length')
  call option(s_coupling,           's.coupling')
  call option(n_length,             'n.length')
  call option(n_spinorbit_a,        'n.spinorbit.a')
  call option(n_spinorbit_b,        'n.spinorbit.b')
  call option(f_length,             'f.length')
  call option(f_exchange_x,         'f.exchange.x')
  call option(f_exchange_y,         'f.exchange.y')
  call option(f_exchange_z,         'f.exchange.z')

  ! Open the output file
  open(newunit=output, file=filename)



  !--------------------------------------------------------------------------------!
  !                           VARIABLE INITIALIZATION                              !
  !--------------------------------------------------------------------------------!

  ! Allocate and initialize the energy array
  allocate(energy_array(energies))
  if (selfconsistent) then
    call energy_range(energy_array, coupling = s_coupling)
  else
    call energy_range(energy_array)
  end if

  ! Calculate the locations of the interfaces
  interfaces(1) = -(s_length + n_length + f_length)/2
  interfaces(2) = interfaces(1) + s_length
  interfaces(3) = interfaces(2) + n_length
  interfaces(4) = interfaces(3) + f_length

  ! Initialize the material layers
  s = superconductor (energy_array, scattering = scattering, thouless = 1/s_length**2, coupling  = s_coupling)
  n = conductor      (energy_array, scattering = scattering, thouless = 1/n_length**2, &
                      spinorbit = spinorbit_xy(alpha = n_spinorbit_a, beta = n_spinorbit_b))
  f = ferromagnet    (energy_array, scattering = scattering, thouless = 1/f_length**2, &
                      exchange  = [f_exchange_x, f_exchange_y, f_exchange_z])

  ! Connect the material layers
  call connect(s, n, conductance, conductance)
  call connect(n, f, conductance, conductance)



  !--------------------------------------------------------------------------------!
  !                           BOOTSTRAPPING PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Write the initial density of states to file
  call write_dos(output)

  ! Bootstrap the system
  do iteration=1,3
    ! Status information
    call print_init

    ! Update the state of the system
    call f%update
    call n%update

    ! Write the density of states to file
    call write_dos(output)
  end do



  !--------------------------------------------------------------------------------!
  !                            CONVERGENCE PROCEDURE                               !
  !--------------------------------------------------------------------------------!

  ! Iterate until convergence
  iteration = 0
  do while (n%difference>10*n%tolerance .or. f%difference>10*f%tolerance .or. (selfconsistent .and. s%difference>10*s%tolerance))
    iteration = iteration + 1

    ! Status information
    call print_main

    ! Update the state of the system
    call f%update
    call n%update
    if (selfconsistent) then
      call s%update
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
    call s % write_dos(output, interfaces(1), interfaces(2))
    call n % write_dos(output, interfaces(2), interfaces(3))
    call f % write_dos(output, interfaces(3), interfaces(4))
    flush(unit=output)
  end subroutine
end program
