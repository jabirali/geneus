! This script calculates the critical temperature of an S/N/F trilayer by performing
! a binary search for the temperature where the gap vanishes numerically. 
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created 2015-07-21
! Updated 2015-07-23

program critical_snf
  use mod_system
  use mod_hybrid
  use mod_critical
  implicit none

  !--------------------------------------------------------------------------------!
  !                            DECLARATION OF VARIABLES                            !
  !--------------------------------------------------------------------------------!

  ! Declare the materials in the hybrid structure
  type(superconductor)    :: s
  type(conductor)         :: n
  type(ferromagnet)       :: f

  ! Declare the materials used for bootstrapping
  type(superconductor)    :: s_boot
  type(conductor)         :: n_boot
  type(ferromagnet)       :: f_boot

  ! Declare the variables used internally by the program
  real(dp),   allocatable :: energy_array(:)
  integer                 :: bisection
  integer                 :: iteration

  ! Declare the parameters that can be modified at runtime
  integer                 :: bisections           = 12
  integer                 :: bootstraps           = 10
  integer                 :: iterations           = 5
  integer                 :: energies             = 600
  real(dp)                :: scattering           = 0.01_dp
  real(dp)                :: conductance          = 0.30_dp
  real(dp)                :: s_gap                = 1e-4_dp
  real(dp)                :: s_length             = 1.00_dp
  real(dp)                :: s_coupling           = 0.20_dp
  real(dp)                :: n_length             = 0.50_dp
  real(dp)                :: n_spinorbit_a        = 0.00_dp
  real(dp)                :: n_spinorbit_b        = 0.00_dp
  real(dp)                :: f_length             = 0.50_dp
  real(dp)                :: f_exchange_x         = 0.00_dp
  real(dp)                :: f_exchange_y         = 0.00_dp
  real(dp)                :: f_exchange_z         = 0.00_dp
  real(dp)                :: temperature_min      = 0.00_dp
  real(dp)                :: temperature_max      = 1.00_dp



  !--------------------------------------------------------------------------------!
  !                           INPUT/OUTPUT PREPARATIONS                            !
  !--------------------------------------------------------------------------------!

  ! Process command line options
  call option
  call option(bisections,           'bisections')
  call option(bootstraps,           'bootstraps')
  call option(iterations,           'iterations')
  call option(energies,             'energies')
  call option(scattering,           'scattering')
  call option(conductance,          'conductance')
  call option(s_gap,                's.gap')
  call option(s_length,             's.length')
  call option(s_coupling,           's.coupling')
  call option(n_length,             'n.length')
  call option(n_spinorbit_a,        'n.spinorbit.a')
  call option(n_spinorbit_b,        'n.spinorbit.b')
  call option(f_length,             'f.length')
  call option(f_exchange_x,         'f.exchange.x')
  call option(f_exchange_y,         'f.exchange.y')
  call option(f_exchange_z,         'f.exchange.z')
  call option(temperature_min,      'temperature.min')
  call option(temperature_max,      'temperature.max')



  !--------------------------------------------------------------------------------!
  !                           VARIABLE INITIALIZATION                              !
  !--------------------------------------------------------------------------------!

  ! Initialize the energy array
  allocate(energy_array(energies))
  call energy_range(energy_array, coupling = s_coupling)

  ! Initialize the backup materials
  s_boot = superconductor(energy_array, coupling = s_coupling)
  n_boot = conductor(energy_array)
  f_boot = ferromagnet(energy_array)

  ! Initialize the material layers
  s = superconductor (energy_array, scattering = scattering, thouless = 1/s_length**2, &
                      coupling  = s_coupling)
  n = conductor      (energy_array, scattering = scattering, thouless = 1/n_length**2, &
                      spinorbit = spinorbit_xy(alpha = n_spinorbit_a, beta = n_spinorbit_b))
  f = ferromagnet    (energy_array, scattering = scattering, thouless = 1/f_length**2, &
                      exchange  = [f_exchange_x, f_exchange_y, f_exchange_z])

  ! Connect the material layers
  call connect(s, n, conductance, conductance)
  call connect(n, f, conductance, conductance)



  !--------------------------------------------------------------------------------!
  !                              BOOTSTRAP PROCEDURE                               !
  !--------------------------------------------------------------------------------!

  ! Initialize all materials to weakly superconducting states
  call init_all(s, cmplx(s_gap,0.0_dp,kind=dp))

  do iteration = 1,bootstraps
    ! Status information
    call print_init

    ! Update the state of the non-superconducting elements
    call f % update
    call n % update
  end do

  ! Save the material states to a backup
  call s % save(s_boot)
  call n % save(n_boot)
  call f % save(f_boot)
  s_boot % greenr = s % greenr
  n_boot % greenr = n % greenr
  f_boot % greenr = f % greenr

  !--------------------------------------------------------------------------------!
  !                  BINARY SEARCH FOR THE CRITICAL TEMPERATURE                    !
  !--------------------------------------------------------------------------------!

  ! Set the system temperature to the midpoint of the search space
  call s%set_temperature( (temperature_min + temperature_max)/2.0_dp )

  ! Perform the binary search
  do bisection = 1,bisections
    ! Status information
    iteration = 0
    call print_main

    ! Load the bootstrapped state
    s % gap = s_gap
    call s % load(s_boot)
    call n % load(n_boot)
    call f % load(f_boot)

    ! Update the superconductor
    call s % update

    ! Perform any additional iterations over the structure
    do iteration = 1,iterations
      ! Status information
      call print_main

      ! Update the system
      call n % update
      call f % update
      call n % update
      call s % update
    end do

    ! Check whether the gap has increased or not, and update the search space accordingly
    if (abs(s%get_gap_mean()/s_gap) >= 1.0_dp) then
      temperature_min = s%get_temperature()
    else
      temperature_max = s%get_temperature()
    end if

    ! Set the system temperature to the midpoint of the updated search space
    call s%set_temperature( (temperature_max + temperature_min)/2.0_dp )
  end do

  ! Print the final results
  call print_final

  !--------------------------------------------------------------------------------!
  !                              CLEANUP PROCEDURE                                 !
  !--------------------------------------------------------------------------------!

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
    write(*,'(a)') '│           BOOTSTRAPPING           │'
    write(*,'(a)') '├───────────────────────────────────┤'
    write(*,'(a,6x,a,i2.2,a,i2.2,7x,a)')                &
      '│','Iteration:     ',iteration,' / ',bootstraps,'│'
    write(*,'(a,6x,a,f8.6,6x,a)')                       &
      '│','Temperature:   ',s%get_temperature(),'│'
    write(*,'(a,6x,a,i2.2,a,i2.2,a,i2.2,6x,a)')         &
      '│','Elapsed time:  ',                            &
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
    write(*,'(a)') '│           BINARY SEARCH           │'
    write(*,'(a)') '├───────────────────────────────────┤'
    write(*,'(a,6x,a,i2.2,a,i2.2,7x,a)')                &
      '│','Bisection:     ',bisection,' / ',bisections,'│'
    write(*,'(a,6x,a,i2.2,a,i2.2,7x,a)')                &
      '│','Iteration:     ',iteration,' / ',iterations,'│'
    write(*,'(a,6x,a,f8.6,6x,a)')                       &
      '│','Temperature:   ',s%get_temperature(),'│'
    write(*,'(a,6x,a,i2.2,a,i2.2,a,i2.2,6x,a)')         &
      '│','Elapsed time:  ',                            &
      int(time/3600.0_sp),':',                          &
      int(mod(time,3600.0_sp)/60.0_sp),':',             &
      int(mod(time,60.0_sp)),                          '│'
    write(*,'(a)') '╘═══════════════════════════════════╛'
  end subroutine

  subroutine print_final
    write(*,'(a)') '                                     '
    write(*,'(a)') '╒═══════════════════════════════════╕'
    write(*,'(a)') '│           FINAL RESULTS           │'
    write(*,'(a)') '├───────────────────────────────────┤'
    write(*,'(a)') '│       Critical temperature:       │'
    write(*,'(a,8x,f18.16,9x,a)') '│', s%get_temperature(), '│'
    write(*,'(a)') '╘═══════════════════════════════════╛'
  end subroutine
end program
