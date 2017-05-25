!> Author:   Jabir Ali Ouassou
!> Date:     2017-02-15
!> Category: Programs
!>
!> This program calculates the charge conductivity of an S/X/N superconducting thin-film structure.

program main
  use :: structure_m
  use :: stdio_m
  use :: math_m
  use :: nambu_m
  use :: matrix_m

  !--------------------------------------------------------------------------------!
  !                           INITIALIZATION PROCEDURE                             !
  !--------------------------------------------------------------------------------!

  ! Declare the superconducting structure
  class(conductor),     pointer :: layer
  type(structure)               :: stack
  type(conductor),      target  :: bulk_a
  type(superconductor), target  :: bulk_b

  ! Declare the nonequilibrium matrices
  integer                                      :: n, m, i, j
  real(wp),    pointer,     dimension(:)       :: energy
  real(wp),    pointer,     dimension(:)       :: location
  complex(wp), allocatable, dimension(:,:,:,:) :: diffusion
  complex(wp), allocatable, dimension(:,:,:)   :: boundary_aa
  complex(wp), allocatable, dimension(:,:,:)   :: boundary_ab
  real(wp),    allocatable, dimension(:)       :: voltage
  real(wp),    allocatable, dimension(:)       :: conductivity

  ! Declare program control parameters
  real(wp), parameter :: threshold = 1e-2
  real(wp), parameter :: tolerance = 1e-8

  ! Redefine stdout and stderr 
  stdout = output('output.log')
  stderr = output('error.log')

  !--------------------------------------------------------------------------------!
  !                            CONSTRUCTION PROCEDURE                              !
  !--------------------------------------------------------------------------------!

  ! Construct the superconducting structure from a configuration file
  stack = structure()
  if (stack % materials() > 1) then
    call error('The material stack should only consist of one layer.')
  end if

  ! Make a pointer to the central layer and its properties
  select type(m => stack % a)
    class is (conductor)
      layer    => m
      energy   => m % energy
      location => m % location
    class default
      call error('This program only supports conductor-class materials.')
  end select

  ! Make sure we're using supported boundary conditions
  if (layer % transparent_b .or. layer % reflecting_b .or. layer % secondorder_b /= 0) then
    call warning('Only first-order tunneling boundary conditions is currently supported by this program.')
  end if

  ! Make the left interface of the central layer transparent
  call layer % conf('transparent_a', 'T')

  ! Construct the surrounding bulk material layers
  bulk_a = conductor()
  bulk_b = superconductor()

  ! Connect the bulk materials to the stack
  layer % material_a => bulk_a
  layer % material_b => bulk_b

  ! Disable the bulk materials from updates
  call bulk_a % conf('order','0')
  call bulk_b % conf('order','0')

  ! Set the inelastic scattering parameters
  bulk_a % scattering_inelastic = layer % scattering_inelastic
  bulk_b % scattering_inelastic = layer % scattering_inelastic



  !--------------------------------------------------------------------------------!
  !                           EQUILIBRIUM CALCULATION                              !
  !--------------------------------------------------------------------------------!

  ! Non-selfconsistent bootstrap procedure
  call stack % converge(threshold = threshold, bootstrap = .true.)

  ! Selfconsistent convergence procedure
  call stack % converge(threshold = tolerance, posthook = posthook)

  ! Status information
  call finalize_equilibrium

  !--------------------------------------------------------------------------------!
  !                        NONEQUILIBRIUM CALCULATION                              !
  !--------------------------------------------------------------------------------!

  ! Allocate working memory
  allocate(diffusion(  0:7, 0:7, size(energy), size(location)), &
           boundary_aa(0:7, 0:7, size(energy)),                 &
           boundary_ab(0:7, 0:7, size(energy)),                 &
           voltage(     -size(energy):+size(energy)),           &
           conductivity(-size(energy):+size(energy))            )

  ! Calculate the diffusion matrix
  call stage('Diffusion matrices')
  do n=1,size(energy)
    do m=1,size(location)
      block
        ! Declare local variables
        type(nambu) :: retarded
        type(nambu) :: advanced

        ! Extract the retarded and advanced propagators
        retarded = layer % propagator(n,m) % retarded()
        advanced = layer % propagator(n,m) % advanced()

        ! Calculate the diffusion matrix coefficients
        do j=0,7
          do i=0,7
            diffusion(i,j,n,m) = trace(nambuv(i) * nambuv(j) - nambuv(i) * retarded * nambuv(j) * advanced)/8
          end do
        end do

        ! Invert the diffusion matrix
        diffusion(:,:,n,m) = inverse(diffusion(:,:,n,m))
      end block
    end do
  end do

  ! Calculate the boundary matrices
  call stage('Interface matrices')
  do n=1,size(energy)
    block
      ! Declare local variables
      type(nambu) :: retarded_a, advanced_a, keldysh_a
      type(nambu) :: retarded_b, advanced_b, keldysh_b

      ! Extract the retarded and advanced propagators on this side
      retarded_a = layer  % propagator(n,size(location)) % retarded()
      advanced_a = layer  % propagator(n,size(location)) % advanced()

      ! Extract the retarded and advanced propagators on the other side
      retarded_b = bulk_b % propagator(n,1) % retarded()
      advanced_b = bulk_b % propagator(n,1) % advanced()

      ! Calculate the boundary matrices
      do j=0,7
        do i=0,7
          ! Calculate the keldysh propagator placeholders
          keldysh_a = retarded_a*nambuv(j) - nambuv(j)*advanced_a
          keldysh_b = retarded_b*nambuv(j) - nambuv(j)*advanced_b

          ! Calculate the boundary matrix coefficients
          boundary_aa(i,j,n) = (layer%conductance_b) * trace(nambuv(i) * (keldysh_a*R(advanced_b) - R(retarded_b)*keldysh_a))/16
          boundary_ab(i,j,n) = (layer%conductance_b) * trace(nambuv(i) * (keldysh_b*T(advanced_a) - T(retarded_a)*keldysh_b))/16
        end do
      end do
    end block
  end do

  ! Initialize the array of applied voltages
  voltage(0) = 0
  do n=1,size(energy)
    voltage(+n) = +energy(n)
    voltage(-n) = -energy(n)
  end do

  ! Calculate the differential conductivity
  call stage('Conductivity')
  do n=1,size(energy)
    block
      complex(wp), dimension(0:7,0:7) :: integral
      
      ! Perform the integral
      integral = 0
      do m=1,size(location)-1
        integral = integral + (location(m+1)-location(m)) * (diffusion(:,:,n,m+1) + diffusion(:,:,n,m))/2
      end do
      integral = matmul(boundary_aa(:,:,n), inverse(matmul(integral,boundary_aa(:,:,n)) - identity(8)))

      ! Invert the result of the integral and save it
      conductivity(+n) = re(integral(4,4) + integral(4,0)) / (layer % conductance_b)
      conductivity(-n) = re(integral(4,4) - integral(4,0)) / (layer % conductance_b)
    end block
  end do

  ! Estimate the zero-bias conductance peak
  conductivity(0) = (conductivity(+1) + conductivity(-1))/2

  ! Write the results to a file
  call finalize_nonequilibrium



  !--------------------------------------------------------------------------------!
  !                                 SUBROUTINES                                    !
  !--------------------------------------------------------------------------------!

contains
  impure subroutine posthook
    ! Write results to output files.
    use :: structure_m

    ! Density of states
    call stack % write_density('density.dat')
  end subroutine

  impure subroutine stage(str)
    ! Write status information to standard out.
    character(len=*) :: str
    integer,  save   :: n = 0
    integer          :: m
    
    ! Calculate iterator
    n = n + 1

    ! Calculate colour
    select case(n)
      case (:1)
        m = 31
      case (2)
        m = 33
      case (3:)
        m = 32
    end select

    ! Status information
    write(stdout,'(2x,"[1m#",i0,2x,"[",i0,"m",a,"[0m")') n, m, str
    flush(stdout)
  end subroutine

  impure subroutine finalize_equilibrium
    ! Write out the final results.
    use :: stdio_m

    ! Status information
    call status_head('EQUILIBRIUM')
    call status_body('State difference', stack % difference())
    call status_body('Charge violation', stack % chargeviolation())
    call status_foot
  end subroutine

  impure subroutine finalize_nonequilibrium
    ! Write out the final results.
    use :: stdio_m

    ! Status information
    call status_head('CONDUCTIVITY')
    call status_body('Zero-bias peak',   conductivity(1))
    call status_foot

    ! Write out the conductivity
    call dump('conductivity.dat', [voltage, conductivity], ['Voltage     ', 'Conductivity'])

    ! Close output files
    close(stdout)
    close(stderr)
  end subroutine

  pure function T(U)
    ! Calculates the contents of the spin-active boundary condition commutators:
    !   I_a ~ [G_a, T(G_b)],   I_b ~ [T(G_a), G_b]
    ! This is used for the calculation of the boundary coefficient matrices.
    !
    ! @NOTE: This version includes only transmission terms.

    type(nambu), intent(in) :: U
    type(nambu)             :: T
    type(nambu)             :: M
    real(wp)                :: GMR
    real(wp)                :: GT1

    ! Extract the interface conductances from the material
    associate(P => layer % polarization_b)
      GMR  = P/(1 + sqrt(1-P**2))
      GT1  = (1 - sqrt(1-P**2))/(1 + sqrt(1-P**2))
    end associate

    ! Extract the interface magnetizations from the material
    M = layer % M_b

    ! Calculate the result
    T = U + GMR * (M*U+U*M) + GT1 * M*U*M
  end function

  pure function R(U)
    ! Calculates the contents of the spin-active boundary condition commutators:
    !   I_a ~ [G_a, R(G_b)],   I_b ~ [R(G_a), G_b]
    ! This is used for the calculation of the boundary coefficient matrices.
    !
    ! @NOTE: This version includes both transmission and reflection terms.

    type(nambu), intent(in) :: U
    type(nambu)             :: R
    type(nambu)             :: M
    complex(wp)             :: Gphi

    ! Extract the interface conductances from the material
    Gphi = (0,-1) * (layer % spinmixing_b)

    ! Extract the interface magnetizations from the material
    M = layer % M0_b

    ! Calculate the result
    R = T(U) + Gphi*M
  end function
end program
