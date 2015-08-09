! This module defines a set of subroutines and functions that are useful for working with multilayer hybrid structures.
! The procedures include the subroutine 'connect' and 'transparent' for creating interfaces between two class(material)
! objects, and a set of functions for initializing arguments that will be passed to class(material) constructor methods.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-11
! Updated: 2015-08-08

module mod_hybrid
  use mod_system
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
    real(dp),                 intent(in)    :: conductance_a   ! Tunneling conductance of the interface (relative to the left  bulk conductance)
    real(dp),                 intent(in)    :: conductance_b   ! Tunneling conductance of the interface (relative to the right bulk conductance)
    
    ! Update the internal material pointers
    material_a % material_b => material_b
    material_b % material_a => material_a

    ! Update the interface parameters
    material_a % conductance_b = conductance_a
    material_b % conductance_a = conductance_b
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
  !              PROCEDURES FOR CLASS(CONDUCTOR) CONSTRUCTOR ARGUMENTS             !
  !--------------------------------------------------------------------------------!

  pure subroutine energy_range(array, coupling)
    ! Initializes an array of energies,  which can be passed on to class(material) constructor methods.
    ! If the coupling constant 'coupling' is provided,  the array will include values all the way up to
    ! the Debye cutoff cosh(1/coupling), which is appropriate for selfconsistent calculations.  If not,
    ! it will only include energies up to 1.5Δ, which is sufficient for non-selfconsistent calculations.

    real(dp), intent(out)          :: array(:)
    real(dp), optional, intent(in) :: coupling
    integer                        :: n

    ! Initialize the energy array
    if (present(coupling) .and. size(array) >= 600) then
      ! Positive energies from 0.0 to 1.5
      do n = 1,size(array)-300
        array(n) = (n-1) * (1.5_dp/(size(array)-300))
      end do
      ! Positive energies from 1.5 to 4.5
      do n = 1,200
        array(size(array)-300+n) = 1.5_dp + (n-1) * (3.0_dp/200)
      end do
      ! Positive energies from 1.5 to cutoff
      do n = 1,100
        array(size(array)-100+n) = 4.5_dp + n * (cosh(1.0_dp/coupling)-4.5)/100
      end do
    else
      ! Positive energies from 0.0 to 1.5
      do n = 1,size(array)
        array(n) = (n-1) * (1.5_dp/(size(array)-1))
      end do
    end if
  end subroutine

  pure function exchange_xy(strength, angle) result(field)
    ! This function returns a vector that describes an exchange field in the xy-plane,
    ! where the input arguments describe the exchange field using polar coordinates.
    real(dp), intent(in) :: strength
    real(dp), intent(in) :: angle
    real(dp)             :: field(3)

    field(1) = strength * cos(angle)
    field(2) = strength * sin(angle)
    field(3) = 0.0_dp
  end function

  pure function spinorbit_xy(strength, angle, alpha, beta) result(field)
    ! This function returns an SU(2) vector that describes a Rashba--Dresselhaus coupling
    ! in the xy-plane.  The coupling constants can also be provided in polar coordinates.
    real(dp), intent(in), optional :: strength
    real(dp), intent(in), optional :: angle
    real(dp), intent(in), optional :: alpha
    real(dp), intent(in), optional :: beta
    type(spin)                     :: field(3)

    ! Initialize to zero
    field(:) = spin(0)

    ! Explicit Rashba coefficient
    if (present(alpha)) then
      field(1) = (+alpha)*pauli2
      field(2) = (-alpha)*pauli1
    end if

    ! Explicit Dresselhaus coefficient
    if (present(beta)) then
      field(1) = (+beta)*pauli1
      field(2) = (-beta)*pauli2
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

  subroutine print_status(header, iteration, change, phasediff, temperature)
    ! Prints a status message to stdout including iteration number, elapsed time, and physical parameters.
    character(*),       intent(in) :: header
    integer,  optional, intent(in) :: iteration
    real(dp), optional, intent(in) :: change
    real(dp), optional, intent(in) :: phasediff
    real(dp), optional, intent(in) :: temperature
    real(sp)                       :: time
    character(len=33)              :: string
    character(len=4)               :: bar

    ! Determine how much CPU time has elapsed
    call cpu_time(time)

    ! Copy the header to a string of correct size
    string = header

    ! Print the progress information to standard out
    if (unicode) then
    write(*,'(a)') '                                     '
    write(*,'(a)') '╒═══════════════════════════════════╕'
    write(*,'(a)') '│ '         // string //          ' │'
    write(*,'(a)') '├───────────────────────────────────┤'
    bar = '│'
    else
    write(*,'(a)') '                                     '
    write(*,'(a)') '+-----------------------------------+'
    write(*,'(a)') '| '         // string //          ' |'
    write(*,'(a)') '+-----------------------------------+'
    bar = '|'
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
    if (unicode) then
    write(*,'(a)') '╘═══════════════════════════════════╛'
    else
    write(*,'(a)') '+-----------------------------------+'
    end if
  end subroutine
end module
