!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This submodule contains procedures for solving the kinetic equations in a material.

submodule (material_m) kinetic_m
  use :: bvp_m
contains
  module procedure kinetic_update
    !! This subroutine calculates the nonequilibrium propagators of a material from its kinetic equations.

    ! Numerical solver
    type(bvp_sol)                                      :: solver

    ! Propagators
    type(propagator), pointer, dimension(:)            :: Gp
    type(propagator), pointer                          :: Ga, Gb
    type(propagator), target                           :: G0

    ! State vectors
    real(wp), dimension(0:15,size(this%location))      :: vp, up
    real(wp), dimension(0:15)                          :: va, vb

    ! Jacobian matrices
    real(wp), dimension(0:15,0:15,size(this%location)) :: Jp
    real(wp), dimension(0:7,0:15)                      :: Jpa, Jaa
    real(wp), dimension(0:7,0:15)                      :: Jpb, Jbb

    ! Loop variables
    integer                                            :: n, m

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

      ! Construct state vectors from the previous distribution (bulk)
      do m=1,size(this % location)
        vp(:7,m) = Gp(m) % h
        vp(8:,m) = Gp(m) % dh
      end do

      ! Construct state vectors from the previous distribution (left)
      va(:7) = Ga % h
      va(8:) = Ga % dh

      ! Construct state vectors from the previous distribution (right)
      vb(:7) = Gb % h
      vb(8:) = Gb % dh

      ! Construct the Jacobian of the differential equation (bulk)
      do m=1,size(this % location)
        call kinetic_jacobian(this, Gp(m), Jp(:,:,m), this % location(m))
      end do

      ! Construct the Jacobian of the differential equation (left)
      call kinetic_jacobian_a(this, Gp(lbound(Gp,1)), Ga, Jpa, Jaa)

      ! Construct the Jacobian of the differential equation (right)
      call kinetic_jacobian_b(this, Gp(ubound(Gp,1)), Gb, Jpb, Jbb)

      ! Initialize the numerical solver
      solver = bvp_init(16, 8, this%location, vp, max_num_subintervals=(size(this%location)*this%scaling))

      ! Solve the differential equation
      solver = bvp_solver(solver, equation, residual, dfdy=jacobian_equation, dbcdy=jacobian_residual, method=this%method, &
                          error_control=this%control, tol=this%tolerance, trace=this%information, stop_on_fail=.false.)

      ! Check if the calculation succeeded
      if (solver % info /= 0) then
        call error('Failed to converge! This is usually because of an ill-posed problem.')
      end if

      ! Extract the numerical solution
      call bvp_eval(solver, this % location, up)

      ! Update the convergence monitor
      this % difference = max(this % difference, maxval(abs(up - vp)))

      ! Update the nonequilibrium state
      do m=1,size(this % location)
        Gp(m) % h  = up(:7,m)
        Gp(m) % dh = up(8:,m)
      end do

      ! Deallocate workspace memory
      call bvp_terminate(solver)
    end do
  contains
    pure subroutine equation(z, u, f)
      !! Definition of the differential equation du/dz=f(z,u).
      real(wp),                  intent(in)  :: z
      real(wp), dimension(0:15), intent(in)  :: u
      real(wp), dimension(0:15), intent(out) :: f
      real(wp), dimension(0:15,0:15)         :: J

      ! Calculate the Jacobian matrix
      call jacobian_equation(z, u, J)

      ! Calculate the state derivative
      f = matmul(J, u)
    end subroutine

    pure subroutine residual(ua, ub, ra, rb)
      !! Definition of the residual equations r(u)=0.
      real(wp), dimension(0:15), intent(in)  :: ua, ub
      real(wp), dimension(0:7),  intent(out) :: ra, rb

      ! Calculate residuals from the boundary conditions (left)
      ra = matmul(Jpa, ua) - matmul(Jaa, va)

      ! Calculate residuals from the boundary conditions (right)
      rb = matmul(Jpb, ub) - matmul(Jbb, vb)
    end subroutine

    pure subroutine jacobian_equation(z, u, J)
      !! Jacobian of the differential equation J=∂f/∂u.
      real(wp),                       intent(in)  :: z
      real(wp), dimension(0:15),      intent(in)  :: u
      real(wp), dimension(0:15,0:15), intent(out) :: J

      ! Interpolate the Jacobian at this position
      J = interpolate(this % location, Jp, z)
    end subroutine

    pure subroutine jacobian_residual(ua, ub, Ja, Jb)
      !! Jacobian of the residual equations J=∂r/∂u.
      real(wp), dimension(0:15),     intent(in)  :: ua, ub
      real(wp), dimension(0:7,0:15), intent(out) :: Ja, Jb

      Ja = Jpa
      Jb = Jpb
    end subroutine
  end procedure

  pure subroutine kinetic_jacobian(this, Gp, Jp, z)
    !! This function calculates the 16×16 Jacobian of the kinetic differential equation.
    class(material),                intent(in)  :: this
    type(propagator),               intent(in)  :: Gp
    real(wp), dimension(0:15,0:15), intent(out) :: Jp
    real(wp),                       intent(in)  :: z
    real(wp), dimension(0:7,0:7)                :: M, Q, dM, dQ, R

    ! Construct the dissipation matrices
    M  = Gp % dissipation()
    dM = Gp % dissipation_gradient()

    ! Construct the condensate matrices
    Q  = Gp % condensate()
    dQ = Gp % condensate_gradient()

    ! Construct the selfenergy matrices
    call this % kinetic_equation(Gp, R, z)

    ! Construct the Jacobian matrix
    Jp(:7,:7) =  0
    Jp(:7,8:) =  identity(8)
    Jp(8:,:7) = -matmul(inverse(M), dQ + R)
    Jp(8:,8:) = -matmul(inverse(M), dM + Q)
  end subroutine

  pure subroutine kinetic_jacobian_a(this, Gp, Ga, Jpa, Jaa)
    !! This function calculates the 8×16 Jacobian of the kinetic boundary condition (left).
    class(material),               intent(in)  :: this
    type(propagator),              intent(in)  :: Gp,  Ga
    real(wp), dimension(0:7,0:15), intent(out) :: Jpa, Jaa
    real(wp), dimension(0:7,0:7)               :: Mp, Qp, Cp, Ca

    if (this % transparent_a) then
      ! Construct a diagonal Jacobian (this side)
      Jpa(:,:7) = identity(8)
      Jpa(:,8:) = 0

      ! Construct a diagonal Jacobian (other side)
      Jaa(:,:7) = identity(8)
      Jaa(:,8:) = 0
    else
      ! Construct the interior matrices
      Mp = Gp % dissipation()
      Qp = Gp % condensate()

      ! Construct the boundary matrices
      call this % kinetic_equation_a(Gp, Ga, Cp, Ca)

      ! Construct the Jacobian (this side)
      Jpa(:,:7) = Cp + Qp
      Jpa(:,8:) =    + Mp

      ! Construct the Jacobian (other side)
      Jaa(:,:7) = Ca
      Jaa(:,8:) = 0
    end if
  end subroutine

  pure subroutine kinetic_jacobian_b(this, Gp, Gb, Jpb, Jbb)
    !! This function calculates the 8×16 Jacobian of the kinetic boundary condition (left).
    class(material),               intent(in)  :: this
    type(propagator),              intent(in)  :: Gp,  Gb
    real(wp), dimension(0:7,0:15), intent(out) :: Jpb, Jbb
    real(wp), dimension(0:7,0:7)               :: Mp, Qp, Cp, Cb

    if (this % transparent_b) then
      ! Construct a diagonal Jacobian (this side)
      Jpb(:,:7) = identity(8)
      Jpb(:,8:) = 0

      ! Construct a diagonal Jacobian (other side)
      Jbb(:,:7) = identity(8)
      Jbb(:,8:) = 0
    else
      ! Construct the interior matrices
      Mp = Gp % dissipation()
      Qp = Gp % condensate()

      ! Construct the boundary matrices
      call this % kinetic_equation_b(Gp, Gb, Cp, Cb)

      ! Construct the Jacobian (this side)
      Jpb(:,:7) = Cp + Qp
      Jpb(:,8:) =    + Mp

      ! Construct the Jacobian (other side)
      Jbb(:,:7) = Cb
      Jbb(:,8:) = 0
    end if
  end subroutine
end submodule
