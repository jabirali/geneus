!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This submodule defines procedures for solving the diffusion equations in a material.

submodule (material_m) diffusion_m
  use :: bvp_m
contains
  module procedure diffusion_update
    !! This subroutine calculates the equilibrium propagators of a material from its diffusion equation.

    ! Numerical solver
    type(bvp_sol)                                 :: solver

    ! Propagators
    type(propagator), pointer, dimension(:)       :: Gp
    type(propagator), pointer                     :: Ga, Gb
    type(propagator), target                      :: G0

    ! State vectors
    real(wp), dimension(0:31,size(this%location)) :: up, vp

    ! Parameters
    complex(wp)                                   :: e

    ! Loop variables
    integer                                       :: n, m

    ! Initialize the state vectors to zero
    vp = 0
    up = 0

    ! Instantiate a vacuum propagator
    G0 = propagator()

    ! Loop over energies
    do n=size(this%energy),1,-1
      ! Status information
      if (this%information >= 0) then
        write(stdout,'(6x,a,1x,i4,1x,a,1x,i4,1x,a,f0.5,a1)',advance='no') &
          '[',n,'/',size(this%energy),']  E = ',this%energy(n), achar(13)
        flush(stdout)
      end if

      ! Update the propagator pointers (bulk)
      Gp => this % propagator(n,:)

      ! Update the propagator pointers (left)
      if (associated(this % material_a)) then
        Ga => this % material_a % propagator(n, ubound(this % material_a % propagator,2))
      else
        Ga => G0
      end if

      ! Update the propagator pointers (right)
      if (associated(this % material_b)) then
        Gb => this % material_b % propagator(n, lbound(this % material_b % propagator,2))
      else
        Gb => G0
      end if

      ! Construct vectors from the previous state
      do m=1,size(this % location)
        call pack(Gp(m), vp(:,m))
      end do

      ! Calculate the complex normalized energy
      e = cx(this%energy(n), this%scattering)/this%thouless

      ! If the difference between iterations is small, use the results from the previous
      ! iteration as an initial guess. If not, use the results form the previous energy.
      if (this % difference < 0.05_wp) then
        solver = bvp_init(32, 16, this%location, vp, max_num_subintervals=(size(this%location)*this%scaling))
      else
        solver = bvp_init(32, 16, this%location, up, max_num_subintervals=(size(this%location)*this%scaling))
      end if

      ! Solve the differential equation
      solver = bvp_solver(solver, equation, residual, method=this%method, error_control=this%control, &
                          tol=this%tolerance, trace=this%information, stop_on_fail=.false.)

      ! Check if the calculation succeeded
      if (solver % info /= 0) then
        call error('Failed to converge! This is usually because of an ill-posed problem.')
      end if

      ! Extract the numerical solution
      call bvp_eval(solver, this % location, up)

      ! Update the convergence monitor
      this % difference = max(this % difference, maxval(abs(up - vp)))

      ! Update the equilibrium state
      do m=1,size(this % location)
        call unpack(Gp(m), up(:,m))
      end do

      ! Deallocate workspace memory
      call bvp_terminate(solver)
    end do
  contains
    pure subroutine equation(z, up, fp)
      !! Definition of the differential equation du/dz=f(z,u).
      real(wp),                  intent(in)  :: z
      real(wp), dimension(0:31), intent(in)  :: up
      real(wp), dimension(0:31), intent(out) :: fp
      type(propagator)                       :: Gp

      ! Construct a propagator object
      call unpack(Gp, up)

      ! Calculate the second-derivatives of the Riccati parameters
      call this % diffusion_equation(Gp, e, z)

      ! Pack the results into a state vector
      fp( 0: 7) = Gp % dg
      fp( 8:15) = Gp % dgt
      fp(16:23) = Gp % d2g
      fp(24:31) = Gp % d2gt
    end subroutine

    pure subroutine residual(ua, ub, ra, rb)
      !! Definition of the residual equations r(u)=0.
      real(wp), dimension(0:31), intent(in)  :: ua, ub
      real(wp), dimension(0:15), intent(out) :: ra, rb

      type(spin)       :: r, rt
      type(propagator) :: Gp

      ! Unpack the state at the left edge
      call unpack(Gp, ua)

      ! Calculate residuals from the boundary conditions
      if (this % transparent_a) then
        ! Transparent interface
        r  = Gp % g  - Ga % g
        rt = Gp % gt - Ga % gt
      else
        ! Spin-active interface
        call this % diffusion_equation_a(Gp, Ga, r, rt)
      end if

      ! Construct the residual vector
      ra(:7) = r
      ra(8:) = rt

      ! Unpack the state at the right edge
      call unpack(Gp, ub)

      ! Calculate residuals from the boundary conditions
      if (this % transparent_b) then
        ! Transparent interface
        r  = Gp % g  - Gb % g
        rt = Gp % gt - Gb % gt
      else
        ! Spin-active interface
        call this % diffusion_equation_b(Gp, Gb, r, rt)
      end if

      ! Construct the residual vector
      rb(:7) = r
      rb(8:) = rt
    end subroutine

    pure subroutine pack(G, u)
      !! Defines assignment from a propagator object to a real vector.
      type(propagator),          intent(in)  :: G
      real(wp), dimension(0:31), intent(out) :: u

      ! Pack the propagator matrices
      u( 0: 7) = G % g
      u( 8:15) = G % gt
      u(16:23) = G % dg
      u(24:31) = G % dgt
    end subroutine

    pure subroutine unpack(G, u)
      !! Defines assignment from a real vector to a propagator object.
      type(propagator),          intent(inout) :: G
      real(wp), dimension(0:31), intent(in)    :: u

      ! Unpack the vector elements
      G % g   = u( 0: 7) 
      G % gt  = u( 8:15) 
      G % dg  = u(16:23) 
      G % dgt = u(24:31) 

      ! Update normalization matrices
      G % N  = inverse( pauli0 - G%g  * G%gt )
      G % Nt = inverse( pauli0 - G%gt * G%g  )
    end subroutine
  end procedure
end submodule
