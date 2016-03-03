! This program calculates the zero-energy peak in the density of states for a superconductor/normal-metal bilayer with
! a spin-active interface for various interface polarizations and spin-dependent phase shifts. The calculation is done
! for magnetizations along different directions, in order to ensure results are independent of the magnetization angle.

program params_spinactive
  use mod_hybrid
  implicit none
!
!  ! Declare variables
!  type(superconductor) :: s               ! Superconductor
!  type(conductor)      :: f               ! Normal metal
!  real(wp)             :: d               ! Density of states at zero excitation energy
!  integer              :: u               ! Output unit
!  integer              :: n, nmax = 200   ! Polarization counter
!  integer              :: m, mmax = 200   ! Phaseshift counter
!  integer              :: k, kmax = 3     ! Magnetization axis
!
!  ! Construct the individual materials
!  s = superconductor ([ 0.0_wp ])
!  f = conductor      ([ 0.0_wp ])
!
!  ! Disable internal output
!  s % information = -1
!  f % information = -1
!
!  ! Connect the two materials
!  call connect(s, f, 0.3_wp)
!
!  do k = 1,kmax
!    select case (k)
!      case (1)
!        ! Open the output file
!        open(newunit=u, file='params_spinactive_x.dat')
!
!        ! Set the magnetization vector
!        f % magnetization_a = [1,0,0]
!      case(2)
!        ! Open the output file
!        open(newunit=u, file='params_spinactive_y.dat')
!
!        ! Set the magnetization vector
!        f % magnetization_a = [0,1,0]
!      case(3)
!        ! Open the output file
!        open(newunit=u, file='params_spinactive_z.dat')
!
!        ! Set the magnetization vector
!        f % magnetization_a = [0,0,1]
!    end select
!
!    do m = -mmax,mmax
!      ! Update the phaseshift
!      f % spinmixing_a = (2*m)/real(mmax,kind=wp)
!
!      do n = 1-nmax,nmax-1
!        ! Update the polarization
!        f % polarization_a = n/real(nmax,kind=wp)
!
!        ! Reinitialize the material
!        call f % init( gap = (1.00_wp,0.01_wp) )
!
!        ! Update the state
!        call f % update
!
!        ! Calculate the zero-energy peak
!        d = f % propagator(1,1) % get_dos()
!
!        ! Write the result to file
!        write(u,'(f10.6,2x)',advance='no') d
!
!        ! Status information
!        write(*,'(1x,a,a,i1,a,i1,a,i4,a,i4,a,i4,a,i4,a,f10.6)',advance='no') &
!          achar(13),' [ Axis: ',k,' / ',kmax,' ] [ Phase: ', m, ' / ', mmax, ' ] [ Polarization: ',n,' / ',nmax,' ]  D(0) =',d
!      end do
!      write(u,*)
!    end do
!
!    ! Close output file
!    close(unit=u)
!  end do
end program
