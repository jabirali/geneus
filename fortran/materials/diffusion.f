!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This submodule defines procedures for solving the diffusion equations in a material.

submodule (material_m) diffusion_m
  use :: bvp_m
contains
  module procedure diffusion_update
    !! This subroutine calculates the equilibrium propagators of a material from the diffusion equation.
    type(bvp_sol)                               :: sol
    real(wp), dimension(32,size(this%location)) :: u, d
    type(propagator)                            :: a, b
    integer                                     :: n, m
    complex(wp)                                 :: e

    ! Loop over energies
    do n=size(this%energy),1,-1
      ! Status information
      if (this%information >= 0) then
        write(stdout,'(6x,a,1x,i4,1x,a,1x,i4,1x,a,f0.5,a1)',advance='no') &
          '[',n,'/',size(this%energy),']  E = ',this%energy(n), achar(13)
        flush(stdout)
      end if

      ! Construct the state vector from the Riccati parameters and their derivatives
      do m=1,size(this%location)
        call pack(this % propagator(n,m), d(:,m))
      end do

      ! If the difference between iterations is small, use the results from the previous
      ! iteration as an initial guess. If not, use the results form the previous energy.
      if (this % difference < 0.05_wp) then
        u = d
      end if

      ! Calculate the complex normalized energy
      e = cx(this%energy(n), this%scattering)/this%thouless

      ! Update the boundary conditions (left)
      a = propagator()
      if (associated(this % material_a)) then
        associate(other => this % material_a % propagator)
          a = other(n,ubound(other,2))
        end associate
      end if

      ! Update the boundary conditions (right)
      b = propagator()
      if (associated(this % material_b)) then
        associate(other => this % material_b % propagator)
          b = other(n,lbound(other,2))
        end associate
      end if

      ! Initialize bvp_solver
      sol = bvp_init(32, 16, this%location, u, max_num_subintervals=(size(this%location)*this%scaling))

      ! Solve the differential equation
      sol = bvp_solver(sol, ode, bc, method=this%method, error_control=this%control, &
                       tol=this%tolerance, trace=this%information, stop_on_fail=.false.)

      ! Check if the calculation succeeded
      if (sol % info /= 0) then
        call error('Failed to converge! This is usually because of an ill-posed problem.')
      end if

      ! Use the results to update the state
      call bvp_eval(sol, this % location, u)
      do m=1,size(this%location)
        call unpack(this % propagator(n,m), u(:,m))
      end do

      ! Update the difference vector
      d = abs(u - d)

      ! Update the maximal difference since last update
      this % difference = max(this % difference, maxval(d))

      ! Clean up after bvp_solver
      call bvp_terminate(sol)
    end do
  contains
    pure subroutine ode(z, u, f)
      !! Definition of the differential equation du/dz=f(z,u).
      real(wp), intent(in)  :: z
      real(wp), intent(in)  :: u(32)
      real(wp), intent(out) :: f(32)
      type(spin)            :: g, gt, dg, dgt
      type(propagator)      :: p

      ! Extract the Riccati parameters
      g   = u( 1: 8)
      gt  = u( 9:16)
      dg  = u(17:24)
      dgt = u(25:32)

      ! Construct a propagator object
      p = propagator(g, gt, dg, dgt)

      ! Calculate the second-derivatives of the Riccati parameters
      call this % diffusion_equation(p, e, z)

      ! Pack the results into a state vector
      f( 1: 8) = p % dg
      f( 9:16) = p % dgt
      f(17:24) = p % d2g
      f(25:32) = p % d2gt
    end subroutine

    pure subroutine bc(ua, ub, bca, bcb)
      !! Definition of the boundary conditions bca=g(ua) and bcb=g(ub).
      real(wp), intent(in)  :: ua(32)
      real(wp), intent(in)  :: ub(32)
      real(wp), intent(out) :: bca(16)
      real(wp), intent(out) :: bcb(16)

      type(propagator)      :: pa, pb
      type(spin)            :: ra, rta
      type(spin)            :: rb, rtb

      ! State at the left end of the material
      call unpack(pa, ua)

      ! State at the right end of the material
      call unpack(pb, ub)

      ! Calculate residuals from the boundary conditions (left)
      if (this % transparent_a) then
        ! Transparent interface
        ra  = pa % g  - a % g
        rta = pa % gt - a % gt
      else
        ! Customized interface
        call this % diffusion_equation_a(pa, a, ra, rta)
      end if

      ! Calculate residuals from the boundary conditions (right)
      if (this % transparent_b) then
        ! Transparent interface
        rb  = pb % g  - b % g
        rtb = pb % gt - b % gt
      else
        ! Customized interface
        call this % diffusion_equation_b(pb, b, rb, rtb)
      end if

      ! Pack the results into state vectors
      bca(1: 8) = ra
      bca(9:16) = rta
      bcb(1: 8) = rb
      bcb(9:16) = rtb
    end subroutine

    pure subroutine pack(a, b)
      !! Defines assignment from a propagator object to a real vector.
      type(propagator),        intent(in)  :: a
      real(wp), dimension(32), intent(out) :: b

      b( 1: 8) = a % g
      b( 9:16) = a % gt
      b(17:24) = a % dg
      b(25:32) = a % dgt
    end subroutine

    pure subroutine unpack(a, b)
      !! Defines assignment from a real vector to a propagator object.
      type(propagator),        intent(inout) :: a
      real(wp), dimension(32), intent(in)    :: b

      ! Unpack the vector elements
      a % g   = b( 1: 8) 
      a % gt  = b( 9:16) 
      a % dg  = b(17:24) 
      a % dgt = b(25:32) 

      ! Update normalization matrices
      a % N  = inverse( pauli0 - a%g  * a%gt )
      a % Nt = inverse( pauli0 - a%gt * a%g  )
    end subroutine
  end procedure
end submodule
