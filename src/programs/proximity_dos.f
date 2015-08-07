! This program calculates the density of states in a multilayer with a single superconducting component,
! which may be treated either self-consistently or not. The density of states as a function of position
! and energy will be written to a given output file.  Note that the output file may be visualized using
! the accompanying Gnuplot script 'plot/dos.plt'. 
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-08-05
! Updated: 2015-08-07

program proximity_dos
  use mod_hybrid
  implicit none

  !--------------------------------------------------------------------------------!
  !                            DECLARATION OF VARIABLES                            !
  !--------------------------------------------------------------------------------!

  ! Declare the materials in the hybrid structure
  type(superconductor)           :: s
  type(ferromagnet), allocatable :: f(:)

  ! Declare the variables used by the program
  integer                        :: output
  real(dp),          allocatable :: energy_array(:)
  real(dp),          allocatable :: connection(:)
  integer                        :: iteration
  integer                        :: n

  ! Declare global parameters that can be modified at runtime
  character(len=32)              :: filename       = 'dos.dat'
  integer                        :: information    = 0
  logical                        :: selfconsistent = .false.
  integer                        :: layers         = 1
  integer                        :: energies       
  integer                        :: points         = 150
  real(dp)                       :: scattering     = 0.01_dp
  real(dp)                       :: conductance    = 0.30_dp



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

  ! Determine how many ferromagnetic layers there are
  call option(layers, 'layers')
  if (layers <= 0) then
    print *
    print *,'Error: there should be minimum one non-superconducting layer in the structure!'
    stop
  end if
  allocate(f(layers))
  allocate(connection(layers+2))

  ! Determine the number of energies to use
  call option(energies, 'energies')
  allocate(energy_array(energies))

  ! Determine the internal position mesh size
  call option(points, 'points')

  ! Determine the inelastic scattering rate
  call option(scattering, 'scattering')

  ! Determine the interface conductance
  call option(conductance, 'conductance')

  ! Flush information to stdout
  flush(unit=stdout)



  !--------------------------------------------------------------------------------!
  !                       PROCESS SUPERCONDUCTING LAYERS                           !
  !--------------------------------------------------------------------------------!

  block
    ! Declare the input variables
    real(dp) :: gap         = 1.00_dp
    real(dp) :: length      = 1.00_dp
    real(dp) :: coupling    = 0.20_dp
    real(dp) :: temperature = 1e-8_dp

    ! Obtain the command line values
    if (selfconsistent) then
      print *
      call option(gap,         's.gap')
      call option(length,      's.length')
      call option(coupling,    's.coupling')
      call option(temperature, 's.temperature')
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
    s = superconductor(energy_array, scattering = scattering, thouless = 1/length**2, points = points, &
                       coupling = coupling, gap = cmplx(gap, 0, kind=dp))

    ! Set the temperature
    if (selfconsistent) then
      call s % set_temperature(temperature)
    end if

    ! Determine the location of the interface
    connection(1) = -length
    connection(2) =  0.0_dp

    ! Set the superconductor change to zero in non-selfconsistent calculations
    if (.not. selfconsistent) then
      s % difference = 0.0_dp
    end if

    ! Set the information level
    s % information = information

    ! Flush information to stdout
    flush(unit=stdout)
  end block



  !--------------------------------------------------------------------------------!
  !                        PROCESS FERROMAGNETIC LAYERS                            !
  !--------------------------------------------------------------------------------!

  do n=1,size(f)
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
      write(ioname, '(a,i0)') 'f', n

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
      f(n) = ferromagnet(energy_array, scattering = scattering, thouless = 1/length**2, points = points, &
                         spinorbit = spinorbit_xy(alpha = spinorbit_a, beta = spinorbit_b),              &
                         exchange  = [exchange_x, exchange_y, exchange_z])

      ! Initialize it to a superconducting state
      call f(n) % init(gap = cmplx(gap,0,kind=dp))

      ! Connect it to the previous material
      if (n > 1) then
        call connect(f(n-1), f(n), conductance, conductance)
      else
        call connect(s,      f(n), conductance, conductance)
      end if

      ! Determine the location of the interface
      connection(n+2) = connection(n+1) + length

      ! Set the information level
      f(n) % information = information

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
    do iteration=1,3
      ! Status information
      call print_status('         INITIALIZATION', iteration = iteration, change = maxval(f%difference))

      ! Update the state of the system (including edges)
      do n=size(f),1,-1
        call f(n) % update
      end do

      ! Update the state of the system (excluding edges)
      if (size(f) > 2) then
        do n=2,size(f)-1
          call f(n) % update
        end do
      end if

      ! Write the density of states to file
      call write_result(output)
    end do
  end if



  !--------------------------------------------------------------------------------!
  !                            CONVERGENCE PROCEDURE                               !
  !--------------------------------------------------------------------------------!

  iteration = 0
  do while (max(maxval(f%difference/f%tolerance),s%difference/s%tolerance) > 10)
    ! Status information
    call print_status('           CONVERGENCE', iteration = iteration, change = max(maxval(f%difference),s%difference))

    ! Update the ferromagnets right-to-left (including edges)
    do n=size(f),1,-1
      call f(n) % update
    end do

    ! Update the superconductor and adjacent ferromagnet (selfconsistent calculations only)
    if (selfconsistent) then
      call s    % update
      call f(1) % update
    end if

    ! Update the ferromagnets left-to-right (excluding edges)
    if (size(f) > 1) then
      do n=2,size(f)-1
        call f(n) % update
      end do
    end if

    ! Write the density of states to file
    call write_result(output)

    ! Update the counter
    iteration = iteration + 1
  end do

  ! Status information
  call print_status('            CONVERGED', change = max(maxval(f%difference),s%difference))


  !--------------------------------------------------------------------------------!
  !                              CLEANUP PROCEDURE                                 !
  !--------------------------------------------------------------------------------!

  ! Close the output file
  close(unit=output)

  ! Deallocate memory
  deallocate(connection)
  deallocate(f)

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
      call s % write_dos(output, connection(1), connection(2))
    end if

    ! Write out the ferromagnet data in all calculations
    do n = 1,size(f)
      call f(n) % write_dos(output, connection(n+1), connection(n+2))
    end do

    ! Flush the changes to file immediately
    flush(unit=output)
  end subroutine
end program
