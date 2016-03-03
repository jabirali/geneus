! This program calculates the density of states in a multilayer structure, which may consist of one or two
! superconducting components that surround an arbitrarily long stack of ferromagnetic layers. The magnetic
! layers will always be treated selfconsistently, while this behaviour is optional for the superconductors.
! The density of states as a function of position and energy will then be written to an output file.  Note
! that the output file may be visualized using the accompanying Gnuplot script 'plot/density.plt'. 
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-08-05
! Updated: 2015-08-10

program density
  use mod_option
  use mod_hybrid
  implicit none

  !--------------------------------------------------------------------------------!
  !                                GLOBAL VARIABLES                                !
  !--------------------------------------------------------------------------------!

  ! Declare the materials in the hybrid structure
  type(superconductor), allocatable :: r(:)            ! Superconducting reservoirs
  type(superconductor), allocatable :: s(:)            ! Superconducting layers
  type(ferromagnet),    allocatable :: f(:)            ! Ferromagnetic layers

  ! Declare global parameters that can be modified at runtime
  character(len=32)                 :: filename        = 'result'
  integer                           :: information     = 0
  logical                           :: selfconsistent  = .false.
  logical                           :: reservoirs      = .false.
  integer                           :: superconductors = 1
  integer                           :: ferromagnets    = 1
  integer                           :: voltages        = 150
  real(wp)                          :: scattering      = 0.01_wp
  real(wp)                          :: phasediff       = 0.00_wp
  real(wp)                          :: phasestep       = 0.15_wp
  real(wp)                          :: temperature     = 0.00_wp
  real(wp)                          :: cutoff          = 30.0_wp

  ! Declare the variables used by the program
  integer                           :: output_density
  integer                           :: output_conduct
  real(wp),             allocatable :: voltage_array(:)
  real(wp),             allocatable :: connection(:)
  integer                           :: n, m



  !--------------------------------------------------------------------------------!
  !                              COMMAND LINE OPTIONS                              !
  !--------------------------------------------------------------------------------!

  ! Process the generic command line options
  call generic_io

  ! Initialize the voltage array
  voltage_array = [ ((1.5_wp*(n-1)/size(voltage_array)), n=1,size(voltage_array)) ]

  ! Construct the left superconducting layer
  call superconductor_io(s, 1)

  ! Construct the ferromagnetic layers
  do m=1,size(f)
    call ferromagnet_io(f, m)
  end do

  ! Construct the right superconducting layer
  if (size(s) > 1) then
    call superconductor_io(s, 2)
  end if



  !--------------------------------------------------------------------------------!
  !                          INITIALIZATION PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Write the initial density of states to file
  call write_density(output_density)

  ! Initialize the ferromagnets
  if (size(f) > 1) then
    n = 0
    do while (maxval(f%difference) > 0.05)
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
      call write_density(output_density)

      ! Update counter
      n = n + 1
    end do
  end if

  !--------------------------------------------------------------------------------!
  !                             BOOTSTRAP PROCEDURE                                !
  !--------------------------------------------------------------------------------!

  if (phasediff > 0) then
    block
      ! Declare block-local variables
      integer     :: steps
      real(wp)    :: phase
      complex(wp) :: gap, gapt

      ! Loop over the number of intermediate phases
      steps = ceiling(phasediff/phasestep)
      do n=0,steps
        ! Calculate the phase
        if (n == 0) then
          phase = 1e-6
        else
          phase = (n*phasediff)/steps
        end if

        ! Calculate the gap
        gap   = exp(((0.0,-0.5)*pi)*phase)
        gapt  = conjg(gap)

        ! Status information
        call print_status('            BOOTSTRAP', iteration = n, change = maxval(f%difference), phasediff = phase)

        ! Reinitialize the superconductors
        call s(1) % init( gap = gap  )
        call s(2) % init( gap = gapt )

        ! Update the ferromagnets right-to-left (including edges)
        do m=size(f),1,-1
          call f(m) % update
        end do

        ! Update the ferromagnets left-to-right (excluding edges)
        if (size(f) > 1) then
          do m=2,size(f)-1
            call f(m) % update
          end do
        end if

        ! Write the density of states to file
        call write_density(output_density)
      end do
    end block
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

    ! Update the ferromagnets left-to-right (excluding edges)
    if (size(f) > 1) then
      do m=2,size(f)-1
        call f(m) % update
      end do
    end if

    ! Update the superconductors (if selfconsistent)
    if (selfconsistent) then
      call s(1) % update
      if (size(s) > 1) then
        call s(2) % update
      end if
    end if

    ! Write the density of states to file
    call write_density(output_density)

    ! Write the differential conductance to file
    call write_conductance(output_conduct)

    ! Update the counter
    n = n + 1

    ! If this is the only material, stop after one iteration
    if ( .not. selfconsistent .and. size(f) <= 1 ) then
      exit
    end if
  end do

  ! Status information
  call print_status('            CONVERGED', change = max(maxval(f%difference),maxval(s%difference)))


  !--------------------------------------------------------------------------------!
  !                              CLEANUP PROCEDURE                                 !
  !--------------------------------------------------------------------------------!

  ! Close the output file
  close(unit=output_density)
  close(unit=output_conduct)

  ! Deallocate memory
  deallocate(voltage_array)
  deallocate(connection)
  deallocate(f)
  deallocate(s)
  if (allocated(r)) then
    deallocate(r)
  end if
contains

  !--------------------------------------------------------------------------------!
  !                              INPUT PROCEDURES                                  !
  !--------------------------------------------------------------------------------!

  subroutine generic_io
    ! This subroutine processes the generic command line options.

    ! Print out the header
    call print_option

    ! Determine the output file
    call option(filename, 'filename')
    open(newunit=output_density, file=(trim(filename) // '.density.dat'))
    open(newunit=output_conduct, file=(trim(filename) // '.conduct.dat'))

    ! Determine the debug level
    call option(information, 'information')
    if (information < -1 .or. information > 2) then
      print *
      print *,'Error: the information level should be in the range [-1,2]!'
      stop
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
      case (3:)
        print *
        print *,'Error: structures with more than two superconductors are not supported!'
        stop
    end select

    ! Determine how many ferromagnetic layers there are
    call option(ferromagnets, 'ferromagnets')
    if (ferromagnets < 0) then
      print *
      print *,'Error: the number of ferromagnetic layers should be zero or a positive number!'
      stop
    else
      allocate(f(ferromagnets))
      allocate(connection(superconductors+ferromagnets+1))
      if (ferromagnets == 0) then
        selfconsistent = .true.
      end if
    end if

    ! Determine whether to perform a selfconsistent calculation
    call option(selfconsistent, 'selfconsistent')
    if (selfconsistent) then
      if (superconductors > 1) then
        reservoirs = .true.
      end if
    else
      cutoff = 0
    end if

    ! Determine whether to use vacuum or reservoirs as the outer boundaries
    if (selfconsistent) then
      call option(reservoirs, 'reservoirs')
      if (reservoirs) then
        allocate(r(superconductors))
      end if
    end if

    call option(voltages, 'voltages')
    if (voltages < 2) then
      print *
      print *,'Error: minimum two voltages required for the calculations!'
      stop
    end if
    allocate(voltage_array(voltages))

    ! Determine the inelastic scattering rate
    call option(scattering, 'scattering')
    if (scattering < 0) then
      print *
      print *,'Error: the scattering parameter should be a positive number!'
      stop
    end if

    ! Determine the temperature
    if (selfconsistent) then
      call option(temperature, 'temperature')
      if (temperature < 0) then
        print *
        print *,'Error: the temperature should be a positive number!'
        stop
      else
        temperature = temperature + 1e-10
      end if
    end if

    ! Determine the phase difference and step
    if (size(s) > 1 .and. ((reservoirs .and. selfconsistent) .or. (.not. selfconsistent))) then
      call option(phasediff, 'phasediff')
      if (phasediff < 0.0_wp) then
        print *
        print *,'Error: the phase difference should be a positive number!'
        stop
      end if
      call option(phasestep, 'phasestep')
      if (phasestep < 0.0_wp) then
        print *
        print *,'Error: the phase steplength should be a positive number!'
        stop
      end if
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
    real(wp)     :: gap        
    real(wp)     :: length     
    real(wp)     :: depairing
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
    logical      :: reflecting_a
    logical      :: reflecting_b

    ! Set the default values
    gap              = 1.00_wp
    length           = 1.00_wp
    magnetization_a  = 0.00_wp
    magnetization_b  = 0.00_wp
    polarization_a   = 0.00_wp
    polarization_b   = 0.00_wp
    conductance_a    = 0.30_wp
    conductance_b    = 0.30_wp
    spinmixing_a     = 0.00_wp
    spinmixing_b     = 0.00_wp
    spinorbit_a      = 0.00_wp
    spinorbit_b      = 0.00_wp
    reflecting_a     = .false.
    reflecting_b     = .false.

    ! Determine the superconductor name
    write(ioname, '(a,i0)') 's', m

    ! Obtain the command line values
    if (m == 1) then
      print *
      print *,'═══════════════════════════════════'
    end if
    if (selfconsistent) then
      call option(reflecting_a,   trim(ioname) // 'l.reflecting')
      call option(conductance_a,   trim(ioname) // 'l.conductance')
      call option(spinmixing_a,    trim(ioname) // 'l.spinmixing')
      call option(polarization_a,  trim(ioname) // 'l.polarization')
      call option(magnetization_a, trim(ioname) // 'l.magnetization')
      print *,'───────────────────────────────────'
    end if
    call option(gap,               trim(ioname) // '.gap')
    if (selfconsistent) then
      call option(length,          trim(ioname) // '.length')
      call option(spinorbit_a,     trim(ioname) // '.rashba')
      call option(spinorbit_b,     trim(ioname) // '.dresselhaus')
      call option(depairing,       trim(ioname) // '.depairing')
      print *,'───────────────────────────────────'
      call option(reflecting_b,   trim(ioname) // 'r.reflecting')
      call option(conductance_b,   trim(ioname) // 'r.conductance')
      call option(spinmixing_b,    trim(ioname) // 'r.spinmixing')
      call option(polarization_b,  trim(ioname) // 'r.polarization')
      call option(magnetization_b, trim(ioname) // 'r.magnetization')
    end if
    if (m == 2) then
      print *,'═══════════════════════════════════'
      print *
      print *
    end if


    ! Construct the superconductor
    s(m) = superconductor(cutoff, scattering = scattering, thouless = 1/length**2, gap = cmplx(gap,0,kind=wp))
    
    ! Construct and connect a superconducting reservoir
    if (reservoirs) then
      if (m == 1) then
        r(m) = superconductor(cutoff, scattering = scattering, gap = exp(((0.0,-0.5)*pi) * phasediff))
        call transparent(r(m), s(m))
      else if (m == 2) then
        r(m) = superconductor(cutoff, scattering = scattering, gap = exp(((0.0,+0.5)*pi) * phasediff))
        call transparent(s(m), r(m))
      end if
    end if

    ! Set the temperature
    s(m) % temperature = temperature

    ! Set the internal fields
    s(m) % depairing   = depairing
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
    s(m) % reflecting_a    = reflecting_a
    s(m) % reflecting_b    = reflecting_b

    ! Connect it to the previous material and determine the interface location
    if (m == 1) then
      connection(1) = 0.0_wp
      connection(2) = length
    else
      connection(size(connection)) = connection(size(connection)-1) + length
      call connect(f(size(f)), s(m))
    end if

    ! For non-selfconsistent calculations, assert that the material has converged
    if (.not. selfconsistent) then
      s(m) % difference = 0.0_wp
    end if

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
    real(wp)         :: gap
    real(wp)         :: length     
    real(wp)         :: depairing
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
    gap              = 1.00_wp
    length           = 1.00_wp
    depairing        = 0.00_wp
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
    call option(gap,             trim(ioname) // '.gap')
    call option(length,          trim(ioname) // '.length')
    call option(spinorbit_a,     trim(ioname) // '.rashba')
    call option(spinorbit_b,     trim(ioname) // '.dresselhaus')
    call option(depairing,       trim(ioname) // '.depairing')
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
    f(m) = ferromagnet(cutoff, scattering = scattering, thouless = 1/length**2, gap = cmplx(gap,0,kind=wp), exchange = exchange)

    ! Set the internal fields
    f(m) % spinorbit = spinorbit_xy(alpha = spinorbit_a, beta = spinorbit_b)
    f(m) % depairing = depairing

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

    ! Determine the location of the interface
    connection(m+2) = connection(m+1) + length

    ! Set the information level
    f(m) % information = information

    ! Flush information to stdout
    flush(unit=stdout)
  end subroutine



  !--------------------------------------------------------------------------------!
  !                              OUTPUT PROCEDURES                                 !
  !--------------------------------------------------------------------------------!

  subroutine write_density(output)
    ! Saves the density of states as a function of position and energy to a file.
    integer, intent(in) :: output

    ! Go to the start of the file
    rewind(unit=output)

    ! Write out the superconductor data (left)
    call s(1) % write_dos(output, connection(1), connection(2))

    ! Write out the ferromagnet data (all)
    do m = 1,size(f)
      call f(m) % write_dos(output, connection(m+1), connection(m+2))
    end do
 
    ! Write out the superconductor data (right)
    if (size(s) > 1) then
      call s(2) % write_dos(output, connection(size(connection)-1), connection(size(connection)))
    end if

    ! Flush the changes to file immediately
    flush(unit=output)
  end subroutine

  subroutine write_conductance(output)
    ! Saves the tunneling current at the first interface as a function of voltage to a file.
    integer,   intent(in) :: output
    real(wp), allocatable :: conductance(:)
    integer               :: n

    if (size(f) > 0) then
      ! Go to the start of the file
      rewind(unit=output)

      ! Calculate the differential conductance
      conductance = differential_conductance(s(1), f(1), voltage_array, s(1)%temperature)

      ! Write the results to file
      do n=size(voltage_array),1,-1
        write(output,*) -voltage_array(n), conductance(n)
      end do
      do n=1,size(voltage_array),+1
        write(output,*) +voltage_array(n), conductance(n)
      end do

      ! Flush the changes to file immediately
      flush(unit=output)

      ! Deallocate the result array
      deallocate(conductance)
    end if
  end subroutine
end program
