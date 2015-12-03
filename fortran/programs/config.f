program config_test
  use mod_config
  use mod_math
  implicit none

  logical  :: selfconsistent = .false.
  integer  :: ferromagnets   = 10
  integer  :: superconductors = 100
  real(wp) :: magnetization(3) = [10,20,30]
  real(wp) :: magnetizations(3) = [10,20,30]
  real(wp) :: test = 0.0_wp
  integer :: iostat = 0
  integer :: unit   = 0
  character(len=132) :: string = ""

  open(newunit=unit, file='config.ini', iostat=iostat, action='read', status='old')
  if (iostat /= 0) then
    call error('failed to open config file')
  end if
  print *,'[global]'
  call config(unit, 'global',         'selfconsistent',  1, 0, selfconsistent)
  call config(unit, 'global',         'superconductors', 1, 0, superconductors)
  call config(unit, 'global',         'ferromagnets',    1, 0, ferromagnets)
  call config(unit, 'global',         'test',            1, 0, test)
  print *
  print *,'[ferromagnet]'
  call config(unit, 'ferromagnet',    'selfconsistent',  1, 0, selfconsistent)
  call config(unit, 'ferromagnet',    'superconductors', 1, 0, superconductors)
  call config(unit, 'ferromagnet',    'ferromagnets',    1, 0, ferromagnets)
  call config(unit, 'ferromagnet',    'test',            1, 0, test)
  call config(unit, 'ferromagnet',    'magnetization',   1, 0, magnetization)
  call config(unit, 'ferromagnet',    'magnetization',   1, +1, magnetization)
  print *,'[superconductor 1]'
  magnetization = [0,0,0]
  call config(unit, 'superconductor', 'magnetization',   1, -1, magnetization)
  magnetization = [0,0,0]
  call config(unit, 'superconductor', 'magnetization',   1,  0, magnetization)
  magnetization = [0,0,0]
  call config(unit, 'superconductor', 'magnetization',   1, +1, magnetization)
  print *
  print *,'[superconductor 2]'
  call config(unit, 'superconductor', 'magnetization',   2, 0, magnetization)
  call config(unit, 'superconductor', 'iterations',      2, 0, ferromagnets)
  print *
  print *,'[superconductor 3]'
  call config(unit, 'superconductor', 'magnetization',   3, 0, magnetizations)
  call config(unit, 'superconductor', 'iterations',      3, 0, ferromagnets)
  print *
  close(unit=unit)
end program  
