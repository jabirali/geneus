! This script calculates the critical temperature of an S/N/F trilayer by performing
! a binary search for the temperature where the gap vanishes numerically. 
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created 2015-07-21
! Updated 2015-07-23

program test_critical
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

  ! Declare the variables used by the program
  real(dp),   allocatable :: energy_array(:)

  ! Declare the parameters that can be modified at runtime
  integer                 :: bisections           = 12
  integer                 :: iterations           = 6
  integer                 :: energies             = 600
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
  real(dp)                :: temperature_min      = 0.00_dp
  real(dp)                :: temperature_max      = 1.00_dp



  !--------------------------------------------------------------------------------!
  !                           INPUT/OUTPUT PREPARATIONS                            !
  !--------------------------------------------------------------------------------!

  ! Process command line options
  call option
  call option(bisections,           'bisections')
  call option(iterations,           'iterations')
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
  call option(temperature_min,      'temperature_min')
  call option(temperature_max,      'temperature_max')



  !--------------------------------------------------------------------------------!
  !                           VARIABLE INITIALIZATION                              !
  !--------------------------------------------------------------------------------!


  ! Initialize the energy array
  allocate(energy_array(energies))
  call energy_range(energy_array, coupling = s_coupling)

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
  !                  BINARY SEARCH FOR THE CRITICAL TEMPERATURE                    !
  !--------------------------------------------------------------------------------!

  ! Perform the binary search for the critical temperature
  call critical_temperature(s, bisections = bisections, iterations = iterations, lower = temperature_min, upper = temperature_max)

  !--------------------------------------------------------------------------------!
  !                              CLEANUP PROCEDURE                                 !
  !--------------------------------------------------------------------------------!

  ! Deallocate memory
  deallocate(energy_array)
end program
