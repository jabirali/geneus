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
  real(dp),          allocatable :: interfaces(:)
  integer                        :: iteration
  integer                        :: n

  ! Declare global parameters that can be modified at runtime
  character(len=64)              :: filename       = 'dos_proximity.dat'
  logical                        :: selfconsistent = .false.
  integer                        :: layers         = 2
  integer                        :: energies       
  real(dp)                       :: scattering     = 0.01_dp
  real(dp)                       :: conductance    = 0.30_dp



  !--------------------------------------------------------------------------------!
  !                          PROCESS COMMAND LINE OPTIONS                          !
  !--------------------------------------------------------------------------------!

  ! Print out the header
  call option

  ! Determine the output file
  call option(filename, 'filename')
  open(newunit=output, file=filename)

  ! Determine whether to perform a selfconsistent calculation
  call option(selfconsistent, 'selfconsistent')
  if (selfconsistent) then
    energies = 600
  else
    energies = 150
  end if

  ! Determine how many ferromagnetic layers there are
  call option(layers, 'layers')
  if (layers > 1) then
    allocate(f(layers-1))
    allocate(interfaces(layers+1))
    interfaces(1) = 0.0_dp
  else
    print *,'Error: there must be at least two layers in the structure!'
    stop
  end if

  ! Determine the number of energies to use
  call option(energies, 'energies')
  allocate(energy_array(energies))

  ! Determine the inelastic scattering rate
  call option(scattering,   'scattering')

  ! Determine the interface conductance
  call option(conductance,  'conductance')

  ! Process the superconducting layers
  block
    ! Declare the input variables
    real(dp) :: length   = 1.00_dp
    real(dp) :: coupling = 0.20_dp

    ! Obtain the command line values
    if (selfconsistent) then
      call option(length,   's.length')
      call option(coupling, 's.coupling')
    end if

    ! Construct a sufficient energy array
    if (selfconsistent .and. energies>200) then
      call energy_range(energy_array, coupling = coupling)
    else
      call energy_range(energy_array)
    end if

    ! Construct the superconductor
    s = superconductor(energy_array, scattering = scattering, thouless = 1/length**2, coupling  = coupling)

    ! Determine the location of the interface
    interfaces(2) = interfaces(1) + length
  end block

  ! Process the ferromagnetic layers
  do n=1,layers-1
    block
      ! Declare the input variables
      character(len=8) :: ioname
      real(dp)         :: length      = 0.50_dp
      real(dp)         :: exchange_x  = 0.00_dp
      real(dp)         :: exchange_y  = 0.00_dp
      real(dp)         :: exchange_z  = 0.00_dp
      real(dp)         :: spinorbit_a = 0.00_dp
      real(dp)         :: spinorbit_b = 0.00_dp

      ! Determine the ferromagnet name
      write(ioname, '(a,i0)') 'f', n

      ! Obtain the command line values
      call option(length,      trim(ioname) // '.length')
      call option(spinorbit_a, trim(ioname) // '.spinorbit.a')
      call option(spinorbit_b, trim(ioname) // '.spinorbit.b')
      call option(exchange_x,  trim(ioname) // '.exchange.x')
      call option(exchange_y,  trim(ioname) // '.exchange.y')
      call option(exchange_z,  trim(ioname) // '.exchange.z')

      ! Construct the ferromagnet
      f(n) = ferromagnet(energy_array, scattering = scattering, thouless = 1/length**2,     &
                                 spinorbit = spinorbit_xy(alpha = spinorbit_a, beta = spinorbit_b), &
                                 exchange  = [exchange_x, exchange_y, exchange_z])

      ! Connect it to the previous material
      if (n > 1) then
        call connect(f(n-1), f(n), conductance, conductance)
      else
        call connect(s,      f(n), conductance, conductance)
      end if

      ! Determine the location of the interface
      interfaces(n+2) = interfaces(n+1) + length
    end block
  end do



  !--------------------------------------------------------------------------------!
  !                           BOOTSTRAPPING PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Write the initial density of states to file
  call write_dos(output)

  ! Bootstrap the ferromagnets
  !if (layers > 2) then
  !  do iteration=1,3
  !    ! Status information
  !    call print_init

  !    ! Update the state of the system
  !    call f%update
  !    call n%update

  !    ! Write the density of states to file
  !    call write_dos(output)
  !  end do
  !end if



  !--------------------------------------------------------------------------------!
  !                            CONVERGENCE PROCEDURE                               !
  !--------------------------------------------------------------------------------!

  ! Iterate until convergence
  iteration = 0
  do while (.not. convergence(s,f))
    iteration = iteration + 1

    ! Status information
    call print_main

    ! Update the ferromagnets right-to-left (including edges)
    do n=layers-1,1,-1
      call f(n) % update
    end do

    ! Update the superconductor and adjacent ferromagnet (selfconsistent calculations only)
    if (selfconsistent) then
      call s    % update
      call f(1) % update
    end if

    ! Update the ferromagnets left-to-right (excluding edges)
    if (layers > 2) then
      do n=1,layers-2
        call f(n) % update
      end do
    end if

    ! Write the density of states to file
    call write_dos(output)
  end do



  !--------------------------------------------------------------------------------!
  !                              CLEANUP PROCEDURE                                 !
  !--------------------------------------------------------------------------------!

  ! Close the output file
  close(unit=output)

  ! Deallocate memory
  deallocate(energy_array)
  deallocate(interfaces)
  deallocate(f)

contains

  !--------------------------------------------------------------------------------!
  !                             LOGICAL PROCEDURES                                 !
  !--------------------------------------------------------------------------------!

  pure function convergence(s, f) result(r)
    ! Returns whether or not the main loop has converged
    class(superconductor), intent(in)  :: s
    class(ferromagnet),    intent(in)  :: f(:)
    logical                            :: r
    integer                            :: n

    r = ((.not. selfconsistent) .or. (s%difference < 10 * s%tolerance))
    do n=1,size(f)
      r = r .and. (f(n)%difference < 10 * f(n)%tolerance)
    end do
  end function

  !--------------------------------------------------------------------------------!
  !                           INPUT/OUTPUT PROCEDURES                              !
  !--------------------------------------------------------------------------------!

  subroutine print_init
    ! Determine how much CPU time has elapsed
    real(sp) :: time
    call cpu_time(time)

    ! Print the progress information to standard out
    write(*,'(a)') '                                     '
    write(*,'(a)') '╒═══════════════════════════════════╕'
    write(*,'(a)') '│          INITIALIZATION           │'
    write(*,'(a)') '├───────────────────────────────────┤'
    write(*,'(a,5x,a,i8,4x,a)')                         &
      '│','Iteration:        ',iteration,'│'
    write(*,'(a,5x,a,i2.2,a,i2.2,a,i2.2,4x,a)')         &
      '│','Elapsed time:     ',                         &
      int(time/3600.0_sp),':',                          &
      int(mod(time,3600.0_sp)/60.0_sp),':',             &
      int(mod(time,60.0_sp)),                          '│'
    write(*,'(a)') '╘═══════════════════════════════════╛'
  end subroutine

  subroutine print_main
    ! Determine how much CPU time has elapsed
    real(sp) :: time
    call cpu_time(time)

    ! Print the progress information to standard out
    write(*,'(a)') '                                     '
    write(*,'(a)') '╒═══════════════════════════════════╕'
    write(*,'(a)') '│           CONVERGENCE             │'
    write(*,'(a)') '├───────────────────────────────────┤'
    write(*,'(a,5x,a,i8,4x,a)')                         &
      '│','Iteration:        ',iteration,'│'
    !write(*,'(a,5x,a,f8.6,4x,a)')                       &
    !  '│','Maximum change:   ',max(n%difference,f%difference),'│'
    write(*,'(a,5x,a,i2.2,a,i2.2,a,i2.2,4x,a)')         &
      '│','Elapsed time:     ',                         &
      int(time/3600.0_sp),':',                          &
      int(mod(time,3600.0_sp)/60.0_sp),':',             &
      int(mod(time,60.0_sp)),                          '│'
    write(*,'(a)') '╘═══════════════════════════════════╛'
  end subroutine

  subroutine write_dos(output)
    integer, intent(in) :: output

    continue
  !  rewind(unit=output)
  !  call s % write_dos(output, interfaces(1), interfaces(2))
  !  call n % write_dos(output, interfaces(2), interfaces(3))
  !  call f % write_dos(output, interfaces(3), interfaces(4))
  !  flush(unit=output)
  end subroutine
end program
