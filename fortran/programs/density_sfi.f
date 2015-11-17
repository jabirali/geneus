! This program calculates the density of states of a superconductor/ferromagnetic-insulator structure,
! where the interface between the layers is fully reflecting, and has strong ferromagnetic properties.
! The density of states, as a function of position and energy, will then be written to an output file.  
! Note that the output file may be visualized with the accompanying Gnuplot script 'plot/density.plt'.
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-11-16
! Updated: 2015-11-17

program density_sfi
  use mod_option
  use mod_hybrid
  implicit none

  !--------------------------------------------------------------------------------!
  !                            DECLARATION OF VARIABLES                            !
  !--------------------------------------------------------------------------------!

  ! Declare the variables used internally by the program
  type(superconductor)  :: s
  real(wp)              :: energies(800)
  integer               :: unit

  ! Declare the parameters that can be modified at runtime
  real(wp)              :: temperature          = 0.00_wp
  real(wp)              :: exponential          = 0.10_wp
  real(wp)              :: conductance          = 1.00_wp
  real(wp)              :: spinmixing           = 0.00_wp
  real(wp)              :: scattering           = 0.05_wp
  real(wp)              :: length               = 3.00_wp

  ! Process command line options
  write(*,*) 'CONFIGURATION:'
  call option(temperature, 'temperature')
  call option(exponential, 'exponential')
  call option(conductance, 'conductance')
  call option(spinmixing,  'spinmixing')
  call option(scattering,  'scattering')
  call option(length,      'length')



  !--------------------------------------------------------------------------------!
  !                           VARIABLE INITIALIZATION                              !
  !--------------------------------------------------------------------------------!

  ! Initialize the energy array
  call energy_range(energies, coupling = s % coupling)

  ! Initialize the superconductor
  s = superconductor(energies)
  s % thouless    = 1/length**2
  s % temperature = temperature
  s % scattering  = scattering

  ! Initialize the interface
  s % reflecting_b    = .true.
  s % spinmixing_b    = spinmixing
  s % magnetization_b = [0.0_wp, 0.0_wp, 1.0_wp]
  s % conductance_b   = 1000

  ! Open output file
  open(newunit=unit, file='density.dat')



  !--------------------------------------------------------------------------------!
  !                              BOOTSTRAP PROCEDURE                               !
  !--------------------------------------------------------------------------------!

  do while (s%conductance_b > conductance)
    ! Status information
    write(*,'(/,1x,a,f7.2,a)') 'BOOTSTRAP [ conductance: ', s%conductance_b, ' ]'

    ! Loop until weak convergence
    do while (s % difference > 0.05)
      ! Update superconductor
      call s % update

      ! Write density of states
      call output
    end do

    ! Update the interface conductance
    s % conductance_b = max(s%conductance_b/(1+exponential), conductance)

    ! Reset the difference counter
    s % difference = 1
  end do



  !--------------------------------------------------------------------------------!
  !                            CONVERGENCE PROCEDURE                               !
  !--------------------------------------------------------------------------------!

  ! Status information
  write(*,'(/,1x,a)') 'CONVERGENCE'

  ! Loop until strong convergence
  do while ((s%difference/s%tolerance) > 10)
    ! Update superconductor
    call s % update

    ! Write density of states
    call output
  end do



  !--------------------------------------------------------------------------------!
  !                              CLEANUP PROCEDURE                                 !
  !--------------------------------------------------------------------------------!

  ! Close the output file
  close(unit=unit)

contains
  impure subroutine output
    integer :: n

    ! Rewind the output file stream
    rewind(unit)

    ! Write out the density of states as a function of energy
    do n = size(s%energy),1,-1
      write(unit,*) -s%energy(n), s%greenr(n,size(s%location))%get_dos()
    end do
    do n = 1,size(s%energy),+1
      write(unit,*) +s%energy(n), s%greenr(n,size(s%location))%get_dos()
    end do

    ! Flush changes to output file
    flush(unit)
  end subroutine
end program
