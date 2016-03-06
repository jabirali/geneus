

program test
  use mod_hybrid
  use mod_material

  type(superconductor) :: s
  type(ferromagnet)    :: f

  s = superconductor(30.0_wp)
  f = ferromagnet(30.0_wp)
  call connect(s,f)

  call f%conf('temperature', '0.10')
  call f%conf('scattering',  '0.05')
  call f%conf('length',      '0.50')

  call s%conf('temperature', '0.10')
  call s%conf('scattering',  '0.05')
  call s%conf('length',      '0.75')
  call s%conf('coupling',    '0.25')

  write(*,*) f%thouless
  write(*,*) f%temperature
  write(*,*) f%scattering

  call f%update
  call s%update
end program
