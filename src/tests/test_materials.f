
program test_materials
  use mod_hybrid
  !use mod_dos
  integer              :: n
  integer              :: u
  real(dp)             :: erg(600)
  type(conductor)      :: c
  type(superconductor) :: s
  type(ferromagnet)    :: f

  call energy_range(erg, coupling = 0.2_dp)

  s = superconductor(erg, thouless = 1/1.0_dp**2, coupling  = 0.2_dp)
  c = conductor     (erg, thouless = 1/0.2_dp**2, spinorbit = spinorbit_xy(2.0_dp,pi/4))
  f = ferromagnet   (erg, thouless = 1/0.2_dp**2, exchange  = exchange_xy(3.0_dp,-pi/4))

  call connect(s, c, 0.33_dp, 0.33_dp)
  call connect(c, f, 0.33_dp, 0.33_dp)

  !do n=1,size(f%location)
  !  f%exchange(1,n) = 0.1_dp*sin(1.57*f%location(n))
  !  f%exchange(2,n) = 0.1_dp*cos(1.57*f%location(n))
  !  f%exchange(3,n) = 0.1_dp
  !end do

  !call calculate_dos(s, iterations=2, unit=1)

  open(newunit=u, file='test_materials.dat') ! Status (Old? New? Replace?) and position (append?) and action (write?)
  call s%write_dos(u, 0.0_dp, 1.0_dp)
  call c%write_dos(u, 1.0_dp, 2.0_dp)
  call f%write_dos(u, 2.0_dp, 3.0_dp)
  close(unit=u)
 
  do n=1,3
    call f%update
         open(newunit=u, file='test_materials.dat') ! Status (Old? New? Replace?) and position (append?) and action (write?)
         call s%write_dos(u, 0.0_dp, 1.0_dp)
         call c%write_dos(u, 1.0_dp, 2.0_dp)
         call f%write_dos(u, 2.0_dp, 3.0_dp)
         close(unit=u)
    call c%update
         open(newunit=u, file='test_materials.dat') ! Status (Old? New? Replace?) and position (append?) and action (write?)
         call s%write_dos(u, 0.0_dp, 1.0_dp)
         call c%write_dos(u, 1.0_dp, 2.0_dp)
         call f%write_dos(u, 2.0_dp, 3.0_dp)
         close(unit=u)
    call s%update
         open(newunit=u, file='test_materials.dat') ! Status (Old? New? Replace?) and position (append?) and action (write?)
         call s%write_dos(u, 0.0_dp, 1.0_dp)
         call c%write_dos(u, 1.0_dp, 2.0_dp)
         call f%write_dos(u, 2.0_dp, 3.0_dp)
         close(unit=u)
    call c%update
         open(newunit=u, file='test_materials.dat') ! Status (Old? New? Replace?) and position (append?) and action (write?)
         call s%write_dos(u, 0.0_dp, 1.0_dp)
         call c%write_dos(u, 1.0_dp, 2.0_dp)
         call f%write_dos(u, 2.0_dp, 3.0_dp)
         close(unit=u)
  end do

  open(newunit=u, file='test_materials.dat') ! Status (Old? New? Replace?) and position (append?) and action (write?)
  call s%write_dos(u, 0.0_dp, 1.0_dp)
  call c%write_dos(u, 1.0_dp, 2.0_dp)
  call f%write_dos(u, 2.0_dp, 3.0_dp)
  close(unit=u)
end program
