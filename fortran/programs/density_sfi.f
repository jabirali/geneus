! This program calculates the density of states of a superconductor/ferromagnetic-insulator structure,
! where the interface between the layers is fully reflecting, and has strong ferromagnetic properties.
! The density of states, as a function of position and energy, will then be written to an output file.  
! Note that the output file may be visualized with the accompanying Gnuplot script 'plot/density.plt'.
!
! Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-11-16
! Updated: 2015-11-18

program density_sfi
  use mod_stdio
  use mod_option
  use mod_hybrid
  implicit none

  !--------------------------------------------------------------------------------!
  !                            DECLARATION OF VARIABLES                            !
  !--------------------------------------------------------------------------------!

  ! Declare the variables used internally by the program
  type(superconductor)  :: s
  integer               :: unit(4)

  ! Declare the parameters that can be modified at runtime
  real(wp)              :: temperature          = 0.00_wp
  real(wp)              :: exponential          = 0.10_wp
  real(wp)              :: conductance          = 1.00_wp
  real(wp)              :: spinmixing           = 0.00_wp
  real(wp)              :: scattering           = 0.01_wp
  real(wp)              :: length               = 1.00_wp
  real(wp)              :: cutoff               = 30.0_wp
  integer               :: information          = 0

  ! Declare iterators used in do-loops
  integer               :: i

  ! Process command line options
  write(*,*) 'CONFIGURATION:'
  call option(information, 'information')
  call option(temperature, 'temperature')
  call option(exponential, 'exponential')
  call option(conductance, 'conductance')
  call option(spinmixing,  'spinmixing')
  call option(scattering,  'scattering')
  call option(length,      'length')



  !--------------------------------------------------------------------------------!
  !                           VARIABLE INITIALIZATION                              !
  !--------------------------------------------------------------------------------!

  ! Initialize the superconductor
  s = superconductor(cutoff)
  s % thouless    = 1/length**2
  s % temperature = temperature
  s % scattering  = scattering
  s % information = information

  ! Initialize the interface
  s % reflecting_b    = .true.
  s % spinmixing_b    = spinmixing
  s % magnetization_b = [0.0_wp, 0.0_wp, 1.0_wp]
  s % conductance_b   = 1000

  ! Open output files
  open(newunit=unit(1), file='dos_l.dat')
  open(newunit=unit(2), file='dos_c.dat')
  open(newunit=unit(3), file='dos_r.dat')
  open(newunit=unit(4), file='gap.dat')



  !--------------------------------------------------------------------------------!
  !                              BOOTSTRAP PROCEDURE                               !
  !--------------------------------------------------------------------------------!

  do while (s%conductance_b > conductance)
    ! Status information
    write(*,'(/,1x,a,f7.2,a)') 'BOOTSTRAP [ conductance: ', s%conductance_b, ' ]'
    flush(stdout)

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
  flush(stdout)

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

  ! Close the output files
  do i=1,4
    close(unit=unit(i))
  end do

contains
  impure subroutine output
    ! Rewind the output file streams
    do i=1,4
      rewind(unit(i))
    end do

    ! Write out the density of states as a function of energy
    do i = size(s%energy),1,-1
      write(unit(1),*) -s%energy(i), s%propagator(i,1)%get_dos()
      write(unit(2),*) -s%energy(i), s%propagator(i,size(s%location)/2)%get_dos()
      write(unit(3),*) -s%energy(i), s%propagator(i,size(s%location))%get_dos()
    end do
    do i = 1,size(s%energy),+1
      write(unit(1),*) +s%energy(i), s%propagator(i,1)%get_dos()
      write(unit(2),*) +s%energy(i), s%propagator(i,size(s%location)/2)%get_dos()
      write(unit(3),*) +s%energy(i), s%propagator(i,size(s%location))%get_dos()
    end do

    ! Write out the superconducting gap as a function of position
    do i = 1,size(s%location)
      write(unit(4),*) s%location(i), abs(s%gap(i))
    end do

    ! Flush changes to output file
    do i=1,4
      flush(unit(i))
    end do
  end subroutine
end program
