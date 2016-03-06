

program test
  use mod_hybrid
  use mod_material

  type(superconductor) :: s
  type(ferromagnet)    :: f

  s = superconductor(30.0_wp)
  f = ferromagnet(30.0_wp)
  call connect(s,f)

  call f%conf('temperature', '0.1')
  call f%conf('scattering',  '0.1')
  call f%conf('length',      '0.5')

  write(*,*) f%thouless
  write(*,*) f%temperature
  write(*,*) f%scattering

  call f%update
end program
