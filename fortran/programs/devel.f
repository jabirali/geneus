! This program does not perform any important operations by itself, but is
! merely used to test new functions and subroutines as they are developed.
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>

program density
  use mod_hybrid
  implicit none

  type(superconductor) :: s1
  type(superconductor) :: s2
  type(ferromagnet)    :: f
  real(wp)             :: phase, x(3)
  real(wp)             :: m(0:3), s(0:3)
  integer              :: n

  x = [3.0_wp, 0.0_wp, 0.0_wp]
  s1 = superconductor(30.0_wp)
  s2 = superconductor(30.0_wp)
  !sj = superconductor(e, coupling=-0.2_wp)
  f  = ferromagnet(30.0_wp, length = 0.5_wp, exchange = x)
  call connect(s1, f , 0.3_wp)
  call connect(f , s2, 0.3_wp)
  f % information = -1
  f % magnetization_a = [1.0_wp,+1.0_wp,0.0_wp]
  f % magnetization_b = [1.0_wp,-1.0_wp,0.0_wp]
  f % spinmixing_a    = 0.5
  f % spinmixing_b    = 0.5
  f % polarization_a  = 0.5
  f % polarization_b  = 0.5

  do n=0,200
    ! Update the phase
    phase = n/200.0
    call s1 % init( gap = exp(((0.0,-0.5)*pi)*phase) )
    call s2 % init( gap = exp(((0.0,+0.5)*pi)*phase) )
    write(*,*) n

    ! Update the system
    call f % update
    !do while (f %difference > 1e-4)
    !  write(*,*) '.'
    !  call f  % update
    !end do

    ! Calculate the density of states
    !write(*,*) f  % density(1,1)

    ! Calculate the mean current
    m(0) = sum(f%current(0,:))/size(f%current(0,:))
    m(1) = sum(f%current(1,:))/size(f%current(1,:))
    m(2) = sum(f%current(2,:))/size(f%current(2,:))
    m(3) = sum(f%current(3,:))/size(f%current(3,:))
    write(*,*) m

    ! Calculate the max deviation 
    s(0) = maxval(abs(f %current(0,:)-m(0)))
    s(1) = maxval(abs(f %current(1,:)-m(1)))
    s(2) = maxval(abs(f %current(2,:)-m(2)))
    s(3) = maxval(abs(f %current(3,:)-m(3)))
    write(*,*) s

    ! Flush output
    flush(stdout)
  end do

end program 
