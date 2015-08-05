! This script calculates the critical temperature of a bulk superconductor, by performing a binary search for the
! temperature where the gap vanishes numerically. The result should be numerically one in the default unit system.
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created 2015-07-21
! Updated 2015-07-23

program test_critical
  use mod_system
  use mod_hybrid
  use mod_critical
  implicit none

  ! Declare variables and initialize their default values
  type(superconductor)  :: s                     ! Superconductor object
  real(dp), allocatable :: domain(:)             ! Discretized energy domain

  integer               :: energies   = 600      ! Number of energies to use in the discretization
  integer               :: bisections = 10       ! Number of bisections (outer loop in the binary search)
  integer               :: iterations = 2        ! Number of iterations (inner loop in the binary search)
  real(dp)              :: coupling   = 0.200_dp ! BCS coupling constant (dimensionless)
  real(dp)              :: lower      = 0.000_dp ! Lower limit for the critical temperature (relative to the critical temperature of a bulk superconductor)
  real(dp)              :: upper      = 1.500_dp ! Upper limit for the critical temperature (relative to the critical temperature of a bulk superconductor)
  real(dp)              :: length     = 25.00_dp ! Superconductor length (relative to the superconducting coherence length)
  real(dp)              :: scattering = 0.001_dp ! Inelastic scattering (imaginary energy contribution)

  ! Process command line options
  call option
  call option(energies,   'energies')
  call option(bisections, 'bisections')
  call option(iterations, 'iterations')
  call option(coupling,   'coupling')
  call option(lower,      'lower')
  call option(upper,      'upper')
  call option(length,     'length')
  call option(scattering, 'scattering')

  ! Initialize the energy array
  allocate(domain(energies))
  call energy_range(domain, coupling = coupling)

  ! Initialize the superconductor
  s = superconductor(domain, coupling = coupling, thouless = 1/length**2, scattering = scattering)

  ! Perform the binary search for the critical temperature
  call critical_temperature(s, bisections = bisections, iterations = iterations, lower = lower, upper = upper)

  ! Deallocate memory
  deallocate(domain)
end program
