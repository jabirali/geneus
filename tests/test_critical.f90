! This script calculates the critical temperature of a bulk superconductor, by performing a binary search for the
! temperature where the gap vanishes numerically. The result should be numerically one in the given unit system.
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created 2015-07-21
! Updated 2015-07-21

program test_critical
  use module_conductor
  use module_superconductor
  implicit none

  type(superconductor) :: s                       ! Superconductor object
  real(dp)             :: erg(600)                ! Discretized energy range
  real(dp)             :: coupling   =  0.200_dp  ! BCS coupling constant
  real(dp)             :: gap        =  0.001_dp  ! Initial value for the superconducting gap (relative to the zero-temperature bulk value)
  real(dp)             :: lower      =  0.000_dp  ! Lower limit for the critical temperature (relative to the bulk value)
  real(dp)             :: upper      =  1.500_dp  ! Upper limit for the critical temperature (relative to the bulk value)
  real(dp)             :: length     = 20.000_dp  ! Length of the superconductor (relative to the correlation length)
  integer              :: iterations = 12         ! Number of iterations of the binary search
  integer              :: n                       ! Loop variable

  ! Initialize the energy array
  call energy_range_positive(erg, 0.2_dp)

  ! Perform the binary search for the critical temperature
  do n=1,12
    ! Initialize a weakly superconducting state at the current temperature
    s = superconductor(erg, cmplx(gap,0.0_dp,kind=dp), 0.2_dp)
    s%temperature = (upper+lower)/2.0_dp
    s%thouless    = 1.0_dp/length**2

    ! Status information
    write(*,'(a)')            '===================================='
    write(*,'(7x,a,i2,a,i2)') 'ITERATION NUMBER: ', n, '/', 12
    write(*,'(a)')            '===================================='

    ! Update the state of the superconductor
    call s%update

    ! Status information
    write(*,*)
    write(*,'(a,f12.7)') ':: Temperature: ', s%temperature
    write(*,'(a,f12.7)') ':: Gap change:  ', abs(s%get_gap_mean())/gap
    write(*,*)

    ! Check whether the mean gap has increased, and act accordingly
    if (abs(s%get_gap_mean()) >= gap) then
      lower = s%temperature
    else
      upper = s%temperature
    end if
  end do

  ! Print out the final result
  print *,'Final estimate for the critical temperature: ', (upper+lower)/2.0_dp
end program
