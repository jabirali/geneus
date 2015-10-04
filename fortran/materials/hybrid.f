! This module defines a set of subroutines and functions that are useful for working with multilayer hybrid structures.
! The procedures include the subroutine 'connect' and 'transparent' for creating interfaces between two class(material)
! objects, and a set of functions for initializing arguments that will be passed to class(material) constructor methods.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-11
! Updated: 2015-10-04

module mod_hybrid
  use mod_stdio
  use mod_math
  use mod_spin
  use mod_material
  use mod_conductor
  use mod_superconductor
  use mod_ferromagnet
  implicit none
contains

  !--------------------------------------------------------------------------------!
  !                PROCEDURES FOR ASSEMBLING MULTILAYER STRUCTURES                 !
  !--------------------------------------------------------------------------------!

  subroutine connect(material_a, material_b, conductance_a, conductance_b)
    ! This subroutine connects two class(material) materials by a tunneling interface, and may
    ! therefore be used to assemble individual material layers to a multilayer hybrid structure.
    class(conductor), target, intent(inout) :: material_a      ! This object represents the left  material
    class(conductor), target, intent(inout) :: material_b      ! This object represents the right material
    real(wp),       optional, intent(in   ) :: conductance_a   ! Tunneling conductance of the interface (relative to the left  bulk conductance)
    real(wp),       optional, intent(in   ) :: conductance_b   ! Tunneling conductance of the interface (relative to the right bulk conductance)
    
    ! Update the internal material pointers
    material_a % material_b => material_b
    material_b % material_a => material_a

    ! Update the interface parameters
    material_a % transparent_b = .false.
    material_b % transparent_a = .false.

    if (present(conductance_a)) then
      material_a % conductance_b = conductance_a
      if (present(conductance_b)) then
        material_b % conductance_a = conductance_b
      else
        material_b % conductance_a = conductance_a
      end if
    end if
  end subroutine

  subroutine transparent(material_a, material_b)
    ! This subroutine connects two class(material) materials by a transparent interface.
    class(conductor), target, intent(inout) :: material_a      ! This object represents the left  material
    class(conductor), target, intent(inout) :: material_b      ! This object represents the right material
    
    ! Update the internal material pointers
    material_a % material_b => material_b
    material_b % material_a => material_a

    ! Update the interface parameters
    material_a % transparent_b = .true.
    material_b % transparent_a = .true.
  end subroutine



  !--------------------------------------------------------------------------------!
  !               PROCEDURES FOR CALCULATING MULTILAYER INTERACTIONS               !
  !--------------------------------------------------------------------------------!

  function differential_conductance(material_a, material_b, voltage, temperature) result(conductance)
    ! Numerically calculates the differential conductance at a tunneling interface.
    class(material), intent(in) :: material_a
    class(material), intent(in) :: material_b
    real(wp),        intent(in) :: temperature
    real(wp),        intent(in) :: voltage(:)
    real(wp),       allocatable :: conductance(:)
    real(wp),       allocatable :: current(:)
    real(wp),       allocatable :: energy(:)
    real(wp),       allocatable :: dos_a(:), ddos_a(:), idos_a(:)
    real(wp),       allocatable :: dos_b(:), ddos_b(:), idos_b(:)
    integer                     :: n, m, err

    ! Allocate memory
    allocate(conductance(size(voltage)))
    allocate(current(size(voltage)))
    allocate(energy(2*size(voltage)))

    allocate(dos_a(size(material_a%energy)))
    allocate(dos_b(size(material_b%energy)))
    allocate(ddos_a(size(material_a%energy)))
    allocate(ddos_b(size(material_b%energy)))
    allocate(idos_a(size(energy)))
    allocate(idos_b(size(energy)))

    ! Initialize the energy array
    call energy_range(energy)

    ! Calculate the density of states at the interface
    do n = 1,size(material_a%energy)
      dos_a(n) = material_a % greenr(n,ubound(material_a%location,1)) % get_dos()
      dos_b(n) = material_b % greenr(n,lbound(material_b%location,1)) % get_dos()
    end do

    ! Create a PCHIP interpolation of the numerical results above
    call dpchez(size(material_a%energy), material_a%energy, dos_a, ddos_a, .false., 0, 0, err)
    call dpchez(size(material_b%energy), material_b%energy, dos_b, ddos_b, .false., 0, 0, err)

    ! Extract the interpolated density of states in the left material
    call dpchfe(size(material_a%energy), material_a%energy, dos_a, ddos_a, 1, .false., size(energy), energy, idos_a, err)

    ! Calculate the current as a function of voltage by numerical integration
    do n = 1,size(voltage)
      ! Extract the voltage-shifted density of states in the right material
      call dpchfe(size(material_b%energy), material_b%energy, dos_b, ddos_b, 1, .false., &
                  size(energy), abs(energy-voltage(n)), idos_b, err)

      ! Calculate the current for this voltage
      current(n) = 0.0_wp
      do m = 1,size(energy)
        current(n) = current(n) + idos_a(m)*idos_b(m)*(fermi(energy(m)-voltage(n))-fermi(energy(m))) &
                                * (maxval(energy)-minval(energy))/(size(energy)-1)
      end do
    end do

    ! Calculate the differential conductance by numerical differentiation (interior points)
    do n = 2,size(voltage)-1
      conductance(n) = (current(n+1)-current(n-1))/(voltage(n+1)-voltage(n-1))
    end do

    ! Calculate the differential conductance by numerical extrapolation (exterior points)
    conductance(lbound(voltage,1)) = 2*conductance(lbound(voltage,1)+1) - conductance(lbound(voltage,1)+2)
    conductance(ubound(voltage,1)) = 2*conductance(ubound(voltage,1)-1) - conductance(ubound(voltage,1)-2)

   ! Deallocate workspace memory
   deallocate(current)
   deallocate(energy)
   deallocate(dos_a)
   deallocate(dos_b)
   deallocate(ddos_a)
   deallocate(ddos_b)
   deallocate(idos_a)
   deallocate(idos_b)
  contains
    pure function fermi(energy)
      ! Evaluates the Fermi function at the given energy and current temperature.
      real(wp), intent(in) :: energy
      real(wp)             :: fermi

      fermi = 1/(exp(energy/(temperature+1e-16))+1)
    end function
  end function



  !--------------------------------------------------------------------------------!
  !              PROCEDURES FOR CLASS(CONDUCTOR) CONSTRUCTOR ARGUMENTS             !
  !--------------------------------------------------------------------------------!

  pure subroutine energy_range(array, coupling)
    ! Initializes an array of energies,  which can be passed on to class(material) constructor methods.
    ! If the coupling constant 'coupling' is provided,  the array will include values all the way up to
    ! the Debye cutoff cosh(1/coupling), which is appropriate for selfconsistent calculations.  If not,
    ! it will only include energies up to 1.5Δ, which is sufficient for non-selfconsistent calculations.

    real(wp), intent(out)          :: array(:)
    real(wp), optional, intent(in) :: coupling
    integer                        :: n

    ! Initialize the energy array
    if (present(coupling) .and. size(array) >= 600) then
      ! Positive energies from 0.0 to 1.5
      do n = 1,size(array)-300
        array(n) = 1e-6_wp + (n-1) * (1.5_wp/(size(array)-300))
      end do
      ! Positive energies from 1.5 to 4.5
      do n = 1,200
        array(size(array)-300+n) = 1e-6_wp + 1.5_wp + (n-1) * (3.0_wp/200)
      end do
      ! Positive energies from 1.5 to cutoff
      do n = 1,100
        array(size(array)-100+n) = 1e-6_wp + 4.5_wp + n * (cosh(1.0_wp/coupling)-4.5)/100
      end do
    else
      ! Positive energies from 0.0 to 1.5
      do n = 1,size(array)
        array(n) = 1e-6 + (n-1) * (1.5_wp/(size(array)-1))
      end do
    end if
  end subroutine

  pure function exchange_xy(strength, angle) result(field)
    ! This function returns a vector that describes an exchange field in the xy-plane,
    ! where the input arguments describe the exchange field using polar coordinates.
    real(wp), intent(in) :: strength
    real(wp), intent(in) :: angle
    real(wp)             :: field(3)

    field(1) = strength * cos(angle)
    field(2) = strength * sin(angle)
    field(3) = 0.0_wp
  end function

  pure function spinorbit_xy(strength, angle, alpha, beta) result(field)
    ! This function returns an SU(2) vector that describes a Rashba--Dresselhaus coupling
    ! in the xy-plane.  The coupling constants can also be provided in polar coordinates.
    real(wp), intent(in), optional :: strength
    real(wp), intent(in), optional :: angle
    real(wp), intent(in), optional :: alpha
    real(wp), intent(in), optional :: beta
    type(spin)                     :: field(3)

    ! Initialize to zero
    field(:) = spin(0)

    ! Explicit Rashba coefficient
    if (present(alpha)) then
      field(1) = field(1) + (+alpha)*pauli2
      field(2) = field(2) + (-alpha)*pauli1
    end if

    ! Explicit Dresselhaus coefficient
    if (present(beta)) then
      field(1) = field(1) + (+beta)*pauli1
      field(2) = field(2) + (-beta)*pauli2
    end if

    ! Polar parametrization
    if (present(strength) .and. present(angle)) then
      field(1) = field(1) + (+strength) * (cos(angle)*pauli1 + sin(angle)*pauli2)
      field(2) = field(2) + (-strength) * (cos(angle)*pauli2 + sin(angle)*pauli1)
    end if
  end function

  !--------------------------------------------------------------------------------!
  !                             INPUT/OUTPUT PROCEDURES                            !
  !--------------------------------------------------------------------------------!

  subroutine print_status(header, bisection, iteration, change, phasediff, temperature)
    ! Prints a status message to stdout including iteration number, elapsed time, and physical parameters.
    character(*),       intent(in) :: header
    integer,  optional, intent(in) :: bisection
    integer,  optional, intent(in) :: iteration
    real(wp), optional, intent(in) :: change
    real(wp), optional, intent(in) :: phasediff
    real(wp), optional, intent(in) :: temperature
    real(sp)                       :: time
    character(len=33)              :: string
    character(len=4)               :: bar

    ! Determine how much CPU time has elapsed
    call cpu_time(time)

    ! Copy the header to a string of correct size
    string = header

    ! Print the progress information to standard out
    write(*,'(a)') '                                     '
    write(*,'(a)') '╒═══════════════════════════════════╕'
    write(*,'(a)') '│ '         // string //          ' │'
    write(*,'(a)') '├───────────────────────────────────┤'
    bar = '│'
    if (present(bisection)) then
      write(*,'(a,3x,a,i8,3x,a)')                       &
        trim(bar),'Bisection:           ', bisection,   trim(bar)
    end if
    if (present(iteration)) then
      write(*,'(a,3x,a,i8,3x,a)')                       &
        trim(bar),'Iteration:           ', iteration,   trim(bar)
    end if
    if (present(phasediff)) then
      write(*,'(a,3x,a,f8.6,3x,a)')                     &
        trim(bar),'Phase difference:    ', phasediff,   trim(bar)
    end if
    if (present(temperature)) then
      write(*,'(a,3x,a,f8.6,3x,a)')                     &
        trim(bar),'Temperature:         ', temperature, trim(bar)
    end if
    if (present(change)) then
      write(*,'(a,3x,a,f8.6,3x,a)')                     &
        trim(bar),'Maximum change:      ', change,      trim(bar)
    end if
    write(*,'(a,3x,a,i2.2,a,i2.2,a,i2.2,3x,a)')         &
      trim(bar),'Elapsed time:        ',                &
      int(time/3600.0_sp),':',                          &
      int(mod(time,3600.0_sp)/60.0_sp),':',             &
      int(mod(time,60.0_sp)),                           trim(bar)
    write(*,'(a)') '╘═══════════════════════════════════╛'

    ! Flush the progress information to standard out
    flush(unit=stdout)
  end subroutine
end module
