! This program calculates the equilibrium propagators in a superconducting multilayer structure, and saves
! some physical observables such as the density of states and superconducting gap to separate output files.
!
! WORK IN PROGRESS:
!   The development of this program was prompted by the maintenance issues with density.f,
!   and this program is supposed to generalize and eventually replace density.f altogether.
!   The ideal is to make it run based on a config file (using mod_config) instead of command
!   line options, and to eventually reduce code duplication between density.f and critical.f.
!
!   Hopefully, the multilayer structure defined here will also replace and supersede the
!   methods in hybrid.f, and make it easier to construct and prototype new driver programs
!   based on superconducting materials. As well as making it easily extensible for new materials.
!
!   The program should write out files dos, dos_1l, dos_1c, dos_1r, ..., conductance_12, ..., gap_1, etc.
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-12-03
! Updated: 2015-12-05


program equilibrium
  use mod_stdio,  only: error
  use mod_config, only: config
  use mod_multilayer
  implicit none

  ! Setting that control the program
  integer :: unit            = 0
  integer :: iostat          = 0
  integer :: energies        = 800
  integer :: debug           = 0

  ! Computational model of the system
  real(wp), allocatable :: e(:)
  type(multilayer), target :: system

  ! Open configuration file
  open(newunit=unit, file='config.ini', action='read', status='old', iostat=iostat)
  if (iostat /= 0) then
    call error('failed to open configuration file ''config.ini''!')
  end if

  ! Process configuration file
  write(*,'(a)') '[equilibrium]'
  call config(unit, 'equilibrium', 'debug',    1, 0, debug)
  write(*,*)

  ! Computation
  call system % config(unit)
  call system % update
  call system % update

  ! Close configuration file
  close(unit)
end program
