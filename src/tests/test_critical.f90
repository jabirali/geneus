! This script calculates the critical temperature of a bulk superconductor, by performing a binary search for the
! temperature where the gap vanishes numerically. The result should be numerically one in the default unit system.
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created 2015-07-21
! Updated 2015-07-23

program test_critical
  use mod_superconductor
  use mod_multilayer
  use mod_critical
  implicit none

  type(superconductor) :: material
  real(dp)             :: energy(600)
  integer              :: n, m

  ! Initialize the energy array
  call energy_range(energy, coupling = 0.200_dp)

  ! Initialize the superconductor
  material = superconductor(energy, coupling = 0.200_dp, thouless = 0.001_dp, scattering = 0.001_dp)

  ! Perform the binary search for the critical temperature
  call critical_temperature(material, bisections = 10, iterations = 2, lower = 0.0_dp, upper = 1.5_dp)
end program
