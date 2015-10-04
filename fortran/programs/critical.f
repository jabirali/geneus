! This program calculates the critical temperature in a multilayer structure, which consists of a single
! superconductor connected to an arbitrarily long stack of ferromagnetic layers. All these materials are
! of course treated selfconsistently. The critical temperature result will be written to standard output.
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-21
! Updated: 2015-08-14

program critical
  use mod_option
  use mod_hybrid
  implicit none

  !--------------------------------------------------------------------------------!
  !                                GLOBAL VARIABLES                                !
  !--------------------------------------------------------------------------------!

  ! Declare the materials in the hybrid structure
  type(superconductor), allocatable :: s(:)                ! Superconducting layer
  type(ferromagnet),    allocatable :: f(:)                ! Ferromagnetic layers
  type(superconductor), allocatable :: sb(:)               ! Superconducting backups
  type(ferromagnet),    allocatable :: fb(:)               ! Ferromagnetic backups

  ! Declare global parameters that can be modified at runtime
  integer                           :: information  = 0
  integer                           :: bisections   = 12
  integer                           :: iterations   = 6
  integer                           :: ferromagnets = 0
  integer                           :: energies     = 800
  integer                           :: points       = 150
  real(wp)                          :: scattering   = 0.01_wp
  real(wp)                          :: coupling     = 0.20_wp
  real(wp)                          :: initgap      = 1e-5_wp
  real(wp)                          :: minimum      = 0.00_wp
  real(wp)                          :: maximum      = 1.00_wp

  ! Declare the variables used by the program
  real(wp),             allocatable :: energy_array(:)
  integer                           :: n, m, k



  !--------------------------------------------------------------------------------!
  !                              COMMAND LINE OPTIONS                              !
  !--------------------------------------------------------------------------------!

  ! Process the generic command line options
  call generic_io

  ! Initialize the energy array
  call energy_range(energy_array, coupling = coupling)

  ! Construct the left superconducting layer
  call superconductor_io(s, 1)

  ! Construct the ferromagnetic layers
  do m=1,size(f)
    call ferromagnet_io(f, m)
  end do

  ! Deallocate the energy array
  deallocate(energy_array)



  !--------------------------------------------------------------------------------!
  !                          INITIALIZATION PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Initialize the materials at zero temperature
  if (size(f) > 0) then
    ! Disable gap updates for the superconductor
    s(1) % coupling = 0.0_wp

    ! Loop until convergence for a fixed gap and zero temperature
    n = 0
    do while (max(maxval(s%difference/s%tolerance),maxval(f%difference/f%tolerance)) > 10)
      ! Status information
      call print_status('         INITIALIZATION', bisection=0, iteration=n, change = maxval(f%difference), temperature = 0.0_wp)

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

      ! Update the state of the superconductor
      call s(1) % update

      ! Update counter
      n = n + 1
    end do

    ! Reenable gap updates
    s(1) % coupling = coupling

    ! Save the current state of the materials
    sb(1) % gap = s(1) % gap
    call s(1) % save(sb(1))
    do n=1,size(f)
      call f(n) % save(fb(n))
    end do
  end if

  !--------------------------------------------------------------------------------!
  !                           BINARY SEARCH PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Set the system temperature to the miwpoint of the search space
  s(1) % temperature = (minimum + maximum)/2.0_wp

  ! Perform the binary search
  do n = 1,bisections
    ! Status information
    call print_status('           CONVERGENCE', bisection=n, iteration=0, temperature = s(1)%temperature)

    ! Load the states from backup
    if (size(f) > 0) then
      s(1) % gap = sb(1) % gap
      call s(1) % load(sb(1))
      do m = 1,size(f)
        call f(m) % load(fb(m))
      end do
    else
      call s(1) % init( gap = cmplx(initgap,0,kind=wp) )
    end if

    ! Update the superconductor once
    call s(1) % update

    ! Perform any additional iterations over the structure
    if (size(f) > 0) then
      do m = 1,iterations
        ! Status information
        call print_status('           CONVERGENCE', bisection=n, iteration=m, temperature = s(1)%temperature)

        ! Update the ferromagnets right-to-left (including edges)
        do k=size(f),1,-1
          call f(k) % update
        end do

        ! Update the ferromagnets left-to-right (excluding edges)
        if (size(f) > 1) then
          do k=2,size(f)-1
            call f(k) % update
          end do
        end if

        ! Update the superconductor
        call s(1) % update
      end do
    end if

    ! Check whether the gap has increased or not, and update the search space accordingly
    if (abs(s(1)%get_gap_mean()/initgap) >= 1.0_wp) then
      minimum = s(1) % temperature
    else
      maximum = s(1) % temperature
    end if

    ! Set the system temperature to the miwpoint of the updated search space
    s(1) % temperature = (maximum + minimum)/2.0_wp
  end do

  ! Print the final results
  call print_status('            CONVERGED', temperature=s(1)%temperature)


  !--------------------------------------------------------------------------------!
  !                              CLEANUP PROCEDURE                                 !
  !--------------------------------------------------------------------------------!

  ! Deallocate memory
  deallocate(s)
  deallocate(f)
  deallocate(sb)
  deallocate(fb)
contains

  !--------------------------------------------------------------------------------!
  !                              INPUT PROCEDURES                                  !
  !--------------------------------------------------------------------------------!

  subroutine generic_io
    ! This subroutine processes the generic command line options.

    ! Print out the header
    call print_option

    ! Determine the debug level
    call option(information, 'information')
    if (information < -1 .or. information > 2) then
      print *
      print *,'Error: the information level should be in the range [-1,2]!'
      stop
    end if

    ! Allocate memory for the superconducting layer
    allocate(s(1))
    allocate(sb(1))

    ! Determine the number of iterations of the binary search
    call option(bisections, 'bisections')
    if (bisections < 1) then
      print *
      print *,'Error: minimum one bisection required for the binary search!'
      stop
    end if

    ! Determine the number of iterations per temperature
    call option(iterations, 'iterations')
    if (iterations < 1) then
      print *
      print *,'Error: minimum one iteration required per temperature!'
      stop
    end if

    ! Allocate memory for the ferromagentic layers
    call option(ferromagnets, 'ferromagnets')
    if (ferromagnets < 0) then
      print *
      print *,'Error: the number of ferromagnetic layers should be zero or a positive number!'
      stop
    else 
      allocate(f(ferromagnets))
      allocate(fb(ferromagnets))
    end if

    ! Determine the number of energies to use
    call option(energies, 'energies')
    if (energies < 600) then
      print *
      print *,'Error: minimum 600 energies required for self-consistent calculations!'
      stop
    end if
    allocate(energy_array(energies))

    ! Determine the internal position mesh size
    call option(points, 'points')
    if (points < 10) then
      print *
      print *,'Error: minimum 10 points required for the calculations!'
      stop
    end if

    ! Determine the inelastic scattering rate
    call option(scattering, 'scattering')
    if (scattering < 0) then
      print *
      print *,'Error: the scattering parameter should be a positive number!'
      stop
    end if

    ! Determine the initial superconducting gap
    call option(initgap, 'initgap')
    if (initgap <= 0) then
      print *
      print *,'Error: the initial superconducting gap should be a positive number!'
      stop
    else if (initgap >= 1) then
      print *
      print *,'Error: the initial superconducting gap should be less than one!'
      stop
    end if

    ! Determine the lower temperature bound
    call option(minimum, 'minimum')
    if (minimum < 0) then
      print *
      print *,'Error: the minimum temperature should be a positive number!'
      stop
    end if

    ! Determine the upper temperature bound
    call option(maximum, 'maximum')
    if (maximum <= minimum) then
      print *
      print *,'Error: the maximum temperature should be higher than the minimum temperature!'
      stop
    end if

    ! Flush information to stdout
    flush(unit=stdout)
  end subroutine



  !--------------------------------------------------------------------------------!
  !                    INTERACTIVE SUPERCONDUCTOR CONSTRUCTOR                      !
  !--------------------------------------------------------------------------------!

  subroutine superconductor_io(s, m)
    ! Constructs a superconductor object based on command line arguments.
    type(superconductor), intent(inout) :: s(:)
    integer,              intent(in   ) :: m

    ! Declare the input/output variables
    character(8) :: ioname
    real(wp)     :: length     
    real(wp)     :: conductance_a
    real(wp)     :: conductance_b
    real(wp)     :: polarization_a
    real(wp)     :: polarization_b
    real(wp)     :: spinmixing_a
    real(wp)     :: spinmixing_b
    real(wp)     :: spinorbit_a
    real(wp)     :: spinorbit_b
    real(wp)     :: magnetization_a(3)
    real(wp)     :: magnetization_b(3)

    ! Set the default values
    length           = 1.00_wp
    spinorbit_a      = 0.00_wp
    spinorbit_b      = 0.00_wp
    magnetization_a  = 0.00_wp
    magnetization_b  = 0.00_wp
    polarization_a   = 0.00_wp
    polarization_b   = 0.00_wp
    conductance_a    = 0.30_wp
    conductance_b    = 0.30_wp

    ! Determine the superconductor name
    write(ioname, '(a,i0)') 's', m

    ! Obtain the command line values
    if (m == 1) then
      print *
      print *,'═══════════════════════════════════'
    end if
    call option(conductance_a,   trim(ioname) // 'l.conductance')
    call option(spinmixing_a,    trim(ioname) // 'l.spinmixing')
    call option(polarization_a,  trim(ioname) // 'l.polarization')
    call option(magnetization_a, trim(ioname) // 'l.magnetization')
    print *,'───────────────────────────────────'
    call option(length,          trim(ioname) // '.length')
    call option(spinorbit_a,     trim(ioname) // '.rashba')
    call option(spinorbit_b,     trim(ioname) // '.dresselhaus')
    print *,'───────────────────────────────────'
    call option(conductance_b,   trim(ioname) // 'r.conductance')
    call option(spinmixing_b,    trim(ioname) // 'r.spinmixing')
    call option(polarization_b,  trim(ioname) // 'r.polarization')
    call option(magnetization_b, trim(ioname) // 'r.magnetization')
    if (size(f) == 0) then
      print *,'═══════════════════════════════════'
      print *
      print *
    end if

    ! Construct the superconductor
    s(m)  = superconductor(energy_array, scattering = scattering, thouless = 1/length**2, &
                          points = points, coupling = coupling, gap = cmplx(initgap,0,kind=wp))

    ! Construct the backup
    sb(m) = superconductor(energy_array, scattering = scattering, thouless = 1/length**2, &
                           points = points, coupling = coupling, gap = cmplx(initgap,0,kind=wp))
    
    ! Set the internal fields
    s(m) % spinorbit   = spinorbit_xy(alpha = spinorbit_a, beta = spinorbit_b)

    ! Set the interface parameters
    s(m) % conductance_a   = conductance_a
    s(m) % conductance_b   = conductance_b
    s(m) % polarization_a  = polarization_a
    s(m) % polarization_b  = polarization_b
    s(m) % spinmixing_a    = spinmixing_a
    s(m) % spinmixing_b    = spinmixing_b
    s(m) % magnetization_a = magnetization_a
    s(m) % magnetization_b = magnetization_b

    ! Set the information level
    s(m) % information = information

    ! Flush information to stdout
    flush(unit=stdout)
  end subroutine



  !--------------------------------------------------------------------------------!
  !                     INTERACTIVE FERROMAGNET CONSTRUCTOR                        !
  !--------------------------------------------------------------------------------!

  subroutine ferromagnet_io(f, m)
    ! Constructs a ferromagnetic object based on command line options.
    type(ferromagnet), intent(inout) :: f(:)
    integer,           intent(in   ) :: m

    ! Declare the input variables
    character(len=8) :: ioname
    real(wp)         :: length     
    real(wp)         :: conductance_a
    real(wp)         :: conductance_b
    real(wp)         :: polarization_a
    real(wp)         :: polarization_b
    real(wp)         :: spinmixing_a
    real(wp)         :: spinmixing_b
    real(wp)         :: spinorbit_a
    real(wp)         :: spinorbit_b
    real(wp)         :: exchange(3)
    real(wp)         :: magnetization_a(3)
    real(wp)         :: magnetization_b(3)

    ! Set the default values
    length           = 1.00_wp
    exchange         = 0.00_wp
    conductance_a    = 0.30_wp
    conductance_b    = 0.30_wp
    magnetization_a  = 0.00_wp
    magnetization_b  = 0.00_wp
    polarization_a   = 0.00_wp
    polarization_b   = 0.00_wp
    spinorbit_a      = 0.00_wp
    spinorbit_b      = 0.00_wp

    ! Determine the ferromagnet name
    write(ioname, '(a,i0)') 'f', m

    ! Obtain the command line values
    print *,'═══════════════════════════════════'
    call option(conductance_a,   trim(ioname) // 'l.conductance')
    call option(spinmixing_a,    trim(ioname) // 'l.spinmixing')
    call option(polarization_a,  trim(ioname) // 'l.polarization')
    call option(magnetization_a, trim(ioname) // 'l.magnetization')
    print *,'───────────────────────────────────'
    call option(length,          trim(ioname) // '.length')
    call option(spinorbit_a,     trim(ioname) // '.rashba')
    call option(spinorbit_b,     trim(ioname) // '.dresselhaus')
    call option(exchange,        trim(ioname) // '.exchange')
    print *,'───────────────────────────────────'
    call option(conductance_b,   trim(ioname) // 'r.conductance')
    call option(spinmixing_b,    trim(ioname) // 'r.spinmixing')
    call option(polarization_b,  trim(ioname) // 'r.polarization')
    call option(magnetization_b, trim(ioname) // 'r.magnetization')
    if (m == size(f)) then
      print *,'═══════════════════════════════════'
      if (size(s) < 2) then
        print *
        print *
      end if
    end if

    ! Construct the ferromagnet
    f(m)  = ferromagnet(energy_array, scattering = scattering, thouless = 1/length**2, &
                        points = points, gap = cmplx(initgap,0,kind=wp), exchange = exchange)

    ! Construct the backup
    fb(m) = ferromagnet(energy_array, scattering = scattering, thouless = 1/length**2, &
                        points = points, gap = cmplx(initgap,0,kind=wp), exchange = exchange)

    ! Set the internal fields
    f(m) % spinorbit = spinorbit_xy(alpha = spinorbit_a, beta = spinorbit_b)

    ! Set the interface parameters
    f(m) % conductance_a   = conductance_a
    f(m) % conductance_b   = conductance_b
    f(m) % polarization_a  = polarization_a
    f(m) % polarization_b  = polarization_b
    f(m) % spinmixing_a    = spinmixing_a
    f(m) % spinmixing_b    = spinmixing_b
    f(m) % magnetization_a = magnetization_a
    f(m) % magnetization_b = magnetization_b

    ! Connect it to the previous material
    if (m == 1) then
      call connect(s(1),   f(m))
    else
      call connect(f(m-1), f(m))
    end if

    ! Set the information level
    f(m) % information = information

    ! Flush information to stdout
    flush(unit=stdout)
  end subroutine
end program
