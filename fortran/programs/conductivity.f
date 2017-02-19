!> Author:   Jabir Ali Ouassou
!> Date:     2017-02-15
!> Category: Programs
!>
!> This program calculates the charge conductivity of an S/X/N superconducting thin-film structure.
!> 
!> @TODO: Check if we need to generalize to spin-active boundary conditions.

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
  class(material),      pointer :: layer
  type(structure)               :: stack
  type(conductor),      target  :: bulk_a
  type(superconductor), target  :: bulk_b

  ! Declare the nonequilibrium matrices
  integer                                      :: n, m, i, j
  real(wp),    pointer,     dimension(:)       :: energy
  real(wp),    pointer,     dimension(:)       :: location
  complex(wp), allocatable, dimension(:,:,:,:) :: diffusion
  complex(wp), allocatable, dimension(:,:,:)   :: boundary_a
  complex(wp), allocatable, dimension(:,:,:)   :: boundary_b
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
  stack = structure('materials.conf')
  if (stack % materials() > 1) then
    call error('The material stack should only consist of one layer.')
  end if

  ! Make a pointer to the central layer and its properties
  layer    => stack % a
  energy   => layer % energy
  location => layer % location

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



  !--------------------------------------------------------------------------------!
  !                           EQUILIBRIUM CALCULATION                              !
  !--------------------------------------------------------------------------------!

  ! Non-selfconsistent bootstrap procedure
  call stack % converge(threshold = threshold, bootstrap = .true.)

  ! Selfconsistent convergence procedure
  call stack % converge(threshold = tolerance, posthook = posthook)

  !--------------------------------------------------------------------------------!
  !                        NONEQUILIBRIUM CALCULATION                              !
  !--------------------------------------------------------------------------------!

  !@TODO: Write status info while looping over energies

  ! Allocate working memory
  allocate(diffusion( 0:7, 0:7, size(energy), size(location)), &
           boundary_a(0:7, 0:7, size(energy)),                 &
           boundary_b(0:7, 0:7, size(energy)),                 &
           conductivity(size(energy))                          )

  ! Calculate the diffusion matrix
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
        do i=0,7
          do j=0,7
            diffusion(i,j,n,m) = trace(nambuv(i) * nambuv(j) - nambuv(i) * retarded * nambuv(j) * advanced)
          end do
        end do

        ! Diffusion should be: 8I, where I identity.
        ! Retarded, advanced should be: ±(diag(1) - antidiag(0.1i))

        ! Normalize and invert the diffusion matrix 
        diffusion(:,:,n,m) = (layer % conductance_b/sqrt(layer % thouless)) * inv(diffusion(:,:,n,m))

        if (n==1 .and. m==size(location)) then
        print *,'retarded re:'
          do i=1,4
            print *,re(retarded % matrix(i,:))
          end do

        print *,'retarded im:'
          do i=1,4
            print *,aimag(retarded % matrix(i,:))
          end do

        print *,'advanced re:'
          do i=1,4
            print *,re(advanced % matrix(i,:))
          end do

        print *,'advanced im:'
          do i=1,4
            print *,aimag(advanced % matrix(i,:))
          end do

        print *,'M re:'
          do i=0,7
            print *,re(diffusion(i,:,n,m))
          end do

        print *,'M im:'
          do i=0,7
            print *,im(diffusion(i,:,n,m))
          end do
        end if
      end block
    end do
  end do

  ! Calculate the boundary matrices
  do n=1,size(energy)
    block
      ! Declare local variables
      type(nambu) :: retarded_a, advanced_a
      type(nambu) :: retarded_b, advanced_b

      ! Extract the retarded and advanced propagators on this side
      retarded_a = layer  % propagator(n,size(location)) % retarded()
      advanced_a = layer  % propagator(n,size(location)) % advanced()

      ! Extract the retarded and advanced propagators on the other side
      retarded_b = bulk_b % propagator(n,1) % retarded()
      advanced_b = bulk_b % propagator(n,1) % advanced()

      ! Calculate the boundary matrix coefficients
      do i=0,7
        do j=0,7
          boundary_a(i,j,n) = trace(nambuv(i) * (retarded_a*nambuv(j) - nambuv(j)*advanced_a) * advanced_b &
                                  - nambuv(i) * retarded_b * (retarded_a*nambuv(j) - nambuv(j)*advanced_a) )

          boundary_b(i,j,n) = trace(nambuv(i) * (retarded_b*nambuv(j) - nambuv(j)*advanced_b) * advanced_a &
                                  - nambuv(i) * retarded_a * (retarded_b*nambuv(j) - nambuv(j)*advanced_b) )
        end do
      end do
    end block
  end do

  ! Calculate the differential conductivity
  do n=1,size(energy)
    block
      complex(wp), dimension(0:7,0:7) :: integral
      
      ! Perform the integral of the inverse diffusion matrix
      integral = inv(boundary_a(:,:,n))
      do m=1,size(location)-1
        integral = integral - (location(m+1)-location(m)) * (diffusion(:,:,n,m+1)+diffusion(:,:,n,m))/4
      end do
      integral = inv(integral)

      ! Invert the result of the integral and save it
      conductivity(n) = re(integral(4,4) - integral(4,0))/16
    end block
  end do

  ! Write the results to a file
  call finalize





  ! @TODO: iterative
  !
  !    Trapezoid integration:
  !     integral[f(z), z=a, z=b] = 0.5·sum( (z(k+1)-z(k))·(f(k+1)-f(k)) )
  !    
  !    Implementation:
  !     h0(0)  = h_alpha
  !     h0(1:) = h(0) + 0.5*cumsum( (z(1:n)-z(0:n-1)) * (Minv(1:n)-Minv(0:n-1)) )
  !    
  !     h(0)  = h0(0)
  !     h(1:) = h0(1:) + i[A, 0.5*cumsum( (z(1:n)-z(0:n-1)) * (h(1:n)-h(0:n-1)) )]



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

  impure subroutine finalize
    ! Write out the final results.
    use :: stdio_m

    ! Status information
    call status_head('COMPLETE')
    call status_body('State difference', stack % difference())
    call status_body('Charge violation', stack % chargeviolation())
    call status_body('Zero-bias result', conductivity(1))
    call status_foot

    ! Write out the conductivity
    call dump('conductivity.dat', [energy, conductivity], ['Voltage     ', 'Conductivity'])

    ! Close output files
    close(stdout)
    close(stderr)
  end subroutine

  impure function inv(a) result(r)
    ! Calculates the inverse of a square matrix.
    complex(dp), dimension(:,:),      intent(in) :: a
    complex(dp), dimension(size(a,1), size(a,1)) :: r

    external    :: zgetrf, zgetri        ! LAPACK routines
    complex(dp) :: work(size(a,1)*64)    ! Working array
    integer     :: ipiv(size(a,1))       ! Pivot indices
    integer     :: info                  ! Error status

    associate(n => size(a,1), m => size(work))
      ! Make a copy
      r = a

      ! LU factorization
      call zgetrf( n, n, r, n, ipiv, info )

      ! Matrix inversion
      call zgetri( n, r, n, ipiv, work, m, info )
    end associate
  end function
end program
