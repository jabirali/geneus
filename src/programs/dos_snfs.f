! This program

! Created: 

program dos_snfs
  use mod_hybrid
  implicit none

  !--------------------------------------------------------------------------------!
  !                            DECLARATION OF VARIABLES                            !
  !--------------------------------------------------------------------------------!

  ! Declare the materials in the hybrid structure
  type(superconductor)    :: s1
  type(conductor)         :: n
  type(ferromagnet)       :: f
  type(superconductor)    :: s2

  ! Declare the parameters used in the program
  integer                 :: output
  real(dp),   allocatable :: energy_array(:)

  complex(dp)             :: gap_left
  complex(dp)             :: gap_middle
  complex(dp)             :: gap_right
  integer :: nn

  real(dp)                :: interfaces(5)

  ! Declare the parameters that can be modified at runtime
  character(len=64)       :: filename             = 'dos_snfs.dat'
  integer                 :: energies             = 150      
  real(dp)                :: scattering           = 0.01_dp
  real(dp)                :: conductance          = 0.30_dp
  real(dp)                :: phase                = 0.00_dp
  real(dp)                :: s_gap                = 1.00_dp
  real(dp)                :: s_coupling           = 0.20_dp
  real(dp)                :: f_exchange_x         = 0.00_dp
  real(dp)                :: f_exchange_y         = 0.00_dp
  real(dp)                :: f_exchange_z         = 0.00_dp
  real(dp)                :: n_spinorbit_a        = 0.00_dp
  real(dp)                :: n_spinorbit_b        = 0.00_dp
  real(dp)                :: s_length             = 1.00_dp
  real(dp)                :: n_length             = 1.00_dp
  real(dp)                :: f_length             = 1.00_dp



  !--------------------------------------------------------------------------------!
  !                            INPUT/OUTPUT SUBROUTINES                            !
  !--------------------------------------------------------------------------------!

  ! Process command line options
  call option
  call option(filename,             'filename')
  call option(energies,             'energies')
  call option(scattering,           'scattering')
  call option(conductance,          'conductance')
  call option(phase,                'phase')
  call option(s_length,             's.length')
  call option(s_coupling,           's.coupling')
  call option(f_length,             'f.length')
  call option(f_exchange_x,         'f.exchange.x')
  call option(f_exchange_y,         'f.exchange.y')
  call option(f_exchange_z,         'f.exchange.z')
  call option(n_length,             'n.length')
  call option(n_spinorbit_a,        'n.spinorbit.a')
  call option(n_spinorbit_b,        'n.spinorbit.b')

  ! Open the output file
  open(newunit=output, file=filename)



  !--------------------------------------------------------------------------------!
  !                            VARIABLE INITIALIZATION                             !
  !--------------------------------------------------------------------------------!

  ! Allocate and initialize the energy array
  allocate(energy_array(energies))
  call energy_range(energy_array)

  ! Calculate the initial gap used in each layer
  gap_left   = exp(+i*phase*pi/2)
  gap_right  = exp(-i*phase*pi/2)

  ! Calculate the locations of the interfaces
  interfaces(1) = 0.0_dp
  interfaces(2) = interfaces(1) + s_length
  interfaces(3) = interfaces(2) + n_length
  interfaces(4) = interfaces(3) + f_length
  interfaces(5) = interfaces(4) + s_length

  ! Initialize the material layers
  s1 = superconductor (energy_array, scattering = scattering, thouless = 1/s_length**2, coupling  = s_coupling)
  s2 = superconductor (energy_array, scattering = scattering, thouless = 1/s_length**2, coupling  = s_coupling)
  n  = conductor      (energy_array, scattering = scattering, thouless = 1/n_length**2, &
                       spinorbit = spinorbit_xy(alpha = n_spinorbit_a, beta = n_spinorbit_b))
  f  = ferromagnet    (energy_array, scattering = scattering, thouless = 1/f_length**2, &
                       exchange  = [f_exchange_x, f_exchange_y, f_exchange_z])

  ! Connect the material layers
  call connect(s1, n, conductance, conductance)
  call connect(n,  f, conductance, conductance)
  call connect(f, s2, conductance, conductance)

  ! BOOTSTRAPPING: Spend some iterations getting there (TODO: Make parameter)
  ! Status information
  write(*,'(/,a,f7.5,a)') 'Setting phase difference to ', 0.0_dp, ' π'
  ! Update the state of the system
  do nn=0,2
    call f%update
    call n%update
  end do
  do nn=1,15
    ! Status information
    write(*,'(/,a,f7.5,a)') 'Setting phase difference to ', abs((phase*nn)/15), ' π'

    ! Reinitialize the superconductors
    call s1%init( exp((+i*pi*phase*nn)/(2*15)) )
    call s2%init( exp((-i*pi*phase*nn)/(2*15)) )

    ! Update the state of the system
    call f%update
    call n%update

    ! Write the density of states to file
    rewind(unit=output)
    call s1 % write_dos(output, interfaces(1), interfaces(2))
    call n  % write_dos(output, interfaces(2), interfaces(3))
    call f  % write_dos(output, interfaces(3), interfaces(4))
    call s2 % write_dos(output, interfaces(4), interfaces(5))
    flush(unit=output)
  end do


  ! Calculate the i
  nn = 0
  do while (n%difference > n%tolerance .or. f%difference > f%tolerance)
    ! Status information
    nn = nn + 1
    write(*,'(/,a,i0)') 'Starting iteration number ', nn

    ! Update the state of the system
    call f%update
    call n%update

    ! Write the density of states to file
    rewind(unit=output)
    call s1 % write_dos(output, interfaces(1), interfaces(2))
    call n  % write_dos(output, interfaces(2), interfaces(3))
    call f  % write_dos(output, interfaces(3), interfaces(4))
    call s2 % write_dos(output, interfaces(4), interfaces(5))
    flush(unit=output)
  end do

  ! Close the output file
  close(unit=output)

  ! Deallocate memory
  deallocate(energy_array)
end program
