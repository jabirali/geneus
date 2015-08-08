! This program calculates the density of states in a multilayer structure, which may consist of one or two
! superconducting components that surround an arbitrarily long stack of ferromagnetic layers. The magnetic
! layers will always be treated selfconsistently, while this behaviour is optional for the superconductors.
! The density of states as a function of position and energy will then be written to an output file.  Note
! that the output file may be visualized using the accompanying Gnuplot script 'plot/density.plt'. 
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-08-05
! Updated: 2015-08-08

program density
  use mod_hybrid
  implicit none

  !--------------------------------------------------------------------------------!
  !                            DECLARATION OF VARIABLES                            !
  !--------------------------------------------------------------------------------!

  ! Declare the materials in the hybrid structure
  type(superconductor), allocatable :: s(:)
  type(ferromagnet),    allocatable :: f(:)

  ! Declare the variables used by the program
  integer                           :: output
  real(dp),             allocatable :: energy_array(:)
  real(dp),             allocatable :: connection(:)
  integer                           :: n, m

  ! Declare global parameters that can be modified at runtime
  character(len=32)                 :: filename        = 'dos.dat'
  integer                           :: information     = 0
  logical                           :: selfconsistent  = .false.
  integer                           :: superconductors = 1
  integer                           :: ferromagnets    = 1
  integer                           :: energies       
  integer                           :: points          = 150
  real(dp)                          :: scattering      = 0.01_dp
  real(dp)                          :: conductance     = 0.30_dp
  real(dp)                          :: phase           = 0.00_dp



  !--------------------------------------------------------------------------------!
  !                    PROCESS GENERIC COMMAND LINE OPTIONS                        !
  !--------------------------------------------------------------------------------!

  ! Print out the header
  call print_option

  ! Determine the output file
  call option(filename, 'filename')
  open(newunit=output, file=filename)

  ! Determine the debug level
  call option(information, 'information')
  print *

  ! Determine whether to perform a selfconsistent calculation
  call option(selfconsistent, 'selfconsistent')
  if (selfconsistent) then
    energies = 600
  else
    energies = 150
  end if

  ! Determine how many superconducting layers there are
  call option(superconductors, 'superconductors')
  select case (superconductors)
    case (:0)
      print *
      print *,'Error: there should be minimum one superconducting layer in the structure!'
      stop
    case (1:2)
      allocate(s(superconductors))
      if (selfconsistent .and. size(s) == 2) then
        print *
        print *,'Error: self-consistent calculations with two superconducting layers are not supported yet!'
        stop
      end if
    case (3:)
      print *
      print *,'Error: structures with more than two superconductors are not supported!'
      stop
  end select

  ! Determine how many ferromagnetic layers there are
  call option(ferromagnets, 'ferromagnets')
  if (ferromagnets <= 0) then
    print *
    print *,'Error: there should be minimum one non-superconducting layer in the structure!'
    stop
  end if
  allocate(f(ferromagnets))
  allocate(connection(superconductors+ferromagnets+1))

  ! Determine the number of energies to use
  call option(energies, 'energies')
  allocate(energy_array(energies))

  ! Determine the internal position mesh size
  call option(points, 'points')

  ! Determine the inelastic scattering rate
  call option(scattering, 'scattering')

  ! Determine the interface conductance
  call option(conductance, 'conductance')

  ! Determine the phase difference
  if (size(s) > 1) then
    call option(phase, 'phase')
  end if

  ! Flush information to stdout
  flush(unit=stdout)



  !--------------------------------------------------------------------------------!
  !                       PROCESS SUPERCONDUCTING LAYERS                           !
  !--------------------------------------------------------------------------------!

  do m=1,superconductors
    block
      ! Declare the input variables
      character(8) :: ioname
      real(dp)     :: gap        
      real(dp)     :: length     
      real(dp)     :: coupling   
      real(dp)     :: temperature

      ! Set the default values
      gap         = 1.00_dp
      length      = 1.00_dp
      coupling    = 0.20_dp
      temperature = 1e-8_dp

      ! Determine the superconductor name
      write(ioname, '(a,i0)') 's', m

      ! Obtain the command line values
      if (selfconsistent) then
        print *
        call option(gap,         trim(ioname) // '.gap')
        call option(length,      trim(ioname) // '.length')
        call option(coupling,    trim(ioname) // '.coupling')
        call option(temperature, trim(ioname) // '.temperature')
      end if

      ! Construct the energy array
      if (selfconsistent) then
        if (energies < 200) then
          print *
          print *,'Error: minimum 200 energies required for self-consistent calculations!'
          stop
        end if
        call energy_range(energy_array, coupling = coupling)
      else
        if (energies < 1) then
          print *
          print *,'Error: minimum one energy required for the calculations!'
          stop
        end if
        call energy_range(energy_array)
      end if

      ! Construct the superconductor
      s(m) = superconductor(energy_array, scattering = scattering, thouless = 1/length**2, points = points, &
                            coupling = coupling, gap = cmplx(gap, 0, kind=dp))

      ! Set the temperature
      if (selfconsistent) then
        call s(m) % set_temperature(temperature)
      end if

      ! Determine the location of the interface
      if (m == 1) then
        connection(1) = 0.0_dp
        connection(2) = length
      else
        connection(size(connection)) = connection(size(connection)-1) + length
      end if

      ! Set the superconductor change to zero in non-selfconsistent calculations
      if (.not. selfconsistent) then
        s(m) % difference = 0.0_dp
      end if

      ! Set the information level
      s(m) % information = information

      ! Flush information to stdout
      flush(unit=stdout)
    end block
  end do



  !--------------------------------------------------------------------------------!
  !                        PROCESS FERROMAGNETIC LAYERS                            !
  !--------------------------------------------------------------------------------!

  do m=1,size(f)
    block
      ! Declare the input variables
      character(len=8) :: ioname
      real(dp)         :: length     
      real(dp)         :: exchange_x 
      real(dp)         :: exchange_y 
      real(dp)         :: exchange_z 
      real(dp)         :: spinorbit_a
      real(dp)         :: spinorbit_b
      real(dp)         :: gap

      ! Set the default values
      gap         = 1.00_dp
      length      = 1.00_dp
      exchange_x  = 0.00_dp
      exchange_y  = 0.00_dp
      exchange_z  = 0.00_dp
      spinorbit_a = 0.00_dp
      spinorbit_b = 0.00_dp

      ! Determine the ferromagnet name
      write(ioname, '(a,i0)') 'f', m

      ! Obtain the command line values
      print *
      call option(gap,         trim(ioname) // '.gap')
      call option(length,      trim(ioname) // '.length')
      call option(spinorbit_a, trim(ioname) // '.spinorbit.a')
      call option(spinorbit_b, trim(ioname) // '.spinorbit.b')
      call option(exchange_x,  trim(ioname) // '.exchange.x')
      call option(exchange_y,  trim(ioname) // '.exchange.y')
      call option(exchange_z,  trim(ioname) // '.exchange.z')

      ! Construct the ferromagnet
      f(m) = ferromagnet(energy_array, scattering = scattering, thouless = 1/length**2, points = points, &
                         spinorbit = spinorbit_xy(alpha = spinorbit_a, beta = spinorbit_b),              &
                         exchange  = [exchange_x, exchange_y, exchange_z])

      ! Initialize it to a superconducting state
      call f(m) % init(gap = cmplx(gap,0,kind=dp))

      ! Connect it to the previous material
      if (m == 1) then
        call connect(s(1),   f(m), conductance, conductance)
      else
        call connect(f(m-1), f(m), conductance, conductance)
      end if

      ! Determine the location of the interface
      connection(m+2) = connection(m+1) + length

      ! Set the information level
      f(m) % information = information

      ! Flush information to stdout
      flush(unit=stdout)
    end block
  end do

  ! Deallocate the energy array
  deallocate(energy_array)



  !--------------------------------------------------------------------------------!
  !                          INITIALIZATION PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Write the initial density of states to file
  call write_result(output)

  ! Initialize the ferromagnets
  if (size(f) > 1) then
    do n=1,3
      ! Status information
      call print_status('         INITIALIZATION', iteration = n, change = maxval(f%difference))

      ! Update the state of the system (including edges)
      do m=size(f),1,-1
        call f(m) % update
      end do

      ! Update the state of the system (excluding edges)
      if (size(f) > 2) then
        do m=2,size(f)-1
          call f(m) % update
        end do
      end if

      ! Write the density of states to file
      call write_result(output)
    end do
  end if



  !--------------------------------------------------------------------------------!
  !                            CONVERGENCE PROCEDURE                               !
  !--------------------------------------------------------------------------------!

  n = 0
  do while (max(maxval(f%difference/f%tolerance),maxval(s%difference/s%tolerance)) > 10)
    ! Status information
    call print_status('           CONVERGENCE', iteration = n, change = max(maxval(f%difference),maxval(s%difference)))

    ! Update the ferromagnets right-to-left (including edges)
    do m=size(f),1,-1
      call f(m) % update
    end do

    ! Update the left superconductor and adjacent ferromagnet (selfconsistent only)
    if (selfconsistent) then
      call s(1) % update
      call f(1) % update
    end if

    ! Update the ferromagnets left-to-right (excluding edges)
    if (size(f) > 1) then
      do m=2,size(f)-1
        call f(m) % update
      end do
    end if

    ! Update the right superconductor and adjacent ferromagnet (selfconsistent only)
    if (selfconsistent .and. size(s) > 1) then
      call f(size(f)) % update
      call s(2)       % update
    end if

    ! Write the density of states to file
    call write_result(output)

    ! Update the counter
    n = n + 1
  end do

  ! Status information
  call print_status('            CONVERGED', change = max(maxval(f%difference),maxval(s%difference)))


  !--------------------------------------------------------------------------------!
  !                              CLEANUP PROCEDURE                                 !
  !--------------------------------------------------------------------------------!

  ! Close the output file
  close(unit=output)

  ! Deallocate memory
  deallocate(connection)
  deallocate(f)
  deallocate(s)

contains

  !--------------------------------------------------------------------------------!
  !                           INPUT/OUTPUT PROCEDURES                              !
  !--------------------------------------------------------------------------------!

  subroutine write_result(output)
    ! Saves the density of states as a function of position and energy to a file.
    integer, intent(in) :: output
    integer             :: n

    ! Go to the start of the file
    rewind(unit=output)

    ! Write out the superconductor data in self-consistent calculations
    if (selfconsistent) then
      call s(1) % write_dos(output, connection(1), connection(2))
    end if

    ! Write out the ferromagnet data in all calculations
    do m = 1,size(f)
      call f(m) % write_dos(output, connection(m+1), connection(m+2))
    end do

    ! Flush the changes to file immediately
    flush(unit=output)
  end subroutine
end program
