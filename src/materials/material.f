! This module defines the data type 'material',  which models the state of a physical material for a discretized range of
! positions and energies. This is an abstract type, meaning that it is not intended to be instantiated on its own, but is
! intended as a base type for physical materials like conductors, superconductors, and ferromagnets. In other words, this
! type defines the essential data structures and program structure, while the derived subtypes will define actual physics.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-29
! Updated: 2015-07-29

module mod_material
  use mod_system
  use mod_spin
  use mod_green
  implicit none

  ! Type declarations
  type, abstract :: material
    ! These parameters determine the basic physical behaviour of a diffusive material
    real(dp)                                  :: thouless        =  1.00_dp         ! Thouless energy of the material (ratio of the diffusion constant to the squared material length)
    real(dp)                                  :: scattering      =  0.01_dp         ! Imaginary energy term (this models inelastic scattering processes and stabilizes the BVP solver)
    real(dp)                                  :: conductance_a   =  0.00_dp         ! Conductance of the left interface  (relative to the bulk conductance of this material)
    real(dp)                                  :: conductance_b   =  0.00_dp         ! Conductance of the right interface (relative to the bulk conductance of this material)

    ! The physical state of the material is modeled as a discretized range of energies, positions, and quasiclassical Green's functions
    real(dp),                     allocatable :: energy(:)                          ! Discretized domain for the energies
    real(dp),                     allocatable :: location(:)                        ! Discretized domain for the positions
    type(green),                  allocatable :: greenr(:,:)                        ! Discretized values for the Green's function (retarded component)
    type(green),                  allocatable :: greenk(:,:)                        ! Discretized values for the Green's function (Keldysh  component)

    ! Hybrid structures are modeled by a double-linked material list, where these two pointers define the neighbours of the current node
    class(material),                  pointer :: material_a      => null()          ! Material connected to this one at the left  interface (default: null pointer, meaning vacuum)
    class(material),                  pointer :: material_b      => null()          ! Material connected to this one at the right interface (default: null pointer, meaning vacuum)

    ! The package bvp_solver is used to handle differential equations, and will be controlled by the following parameters
    integer                                   :: scaling         =  64              ! Maximal allowed increase in the mesh resolution (range: 2^N, N>1)
    integer                                   :: order           =  4               ! Order of the Runge—Kutta method used by the solver (range: 2, 4, 6)
    integer                                   :: control         =  2               ! Error control method (1: defect, 2: global error, 3: 1 then 2, 4: 1 and 2)
    integer                                   :: information     =  0               ! How much information that should be written to standard out (range: [-1,2])
    real(dp)                                  :: tolerance       =  1e-4_dp         ! Error tolerance (determines the maximum allowed defect or global error)

    ! The following variables are used for input/output purposes, and should be modified by class(material) constructors
    character(len=64)                         :: type_string     =  'MATERIAL'      ! The type string should describe the specific class(material) subtype
  contains
    ! These methods define how to update the physical state of the material
    procedure(init),                 deferred :: init                               ! Initializes the Green's functions
    procedure                                 :: update          => material_update ! Calculates  the Green's functions
    procedure(update),               deferred :: update_prehook                     ! Code to execute before calculating the Green's functions
    procedure(update),               deferred :: update_posthook                    ! Code to execute after  calculating the Green's functions

    ! These methods define the physical equations used by the update methods
    procedure(diffusion_equation),   deferred :: diffusion_equation                 ! Diffusion equation that describes the material
    procedure(interface_equation_a), deferred :: interface_equation_a               ! Boundary condition at the left  interface
    procedure(interface_equation_b), deferred :: interface_equation_b               ! Boundary condition at the right interface
  end type

  ! Interface declarations
  abstract interface
    subroutine init(this, gap)
      ! This interface is used for the deferred procedure init.
      import material, dp

      class(material), intent(inout) :: this
      complex(dp),     intent(in   ) :: gap
    end subroutine
  end interface

  abstract interface
    subroutine update(this)
      ! This interface is used for the deferred procedures update_prehook and update_posthook.
      import material

      class(material), intent(inout) :: this
    end subroutine
  end interface

  abstract interface
    subroutine diffusion_equation(this, e, z, g, gt, dg, dgt, d2g, d2gt)
      ! This interface is used for the deferred procedure diffusion_equation.
      import material, spin, dp

      class(material), intent(in   ) :: this
      type(spin),      intent(in   ) :: g, gt, dg, dgt
      type(spin),      intent(inout) :: d2g, d2gt
      complex(dp),     intent(in   ) :: e
      real(dp),        intent(in   ) :: z
    end subroutine
  end interface

  abstract interface
    subroutine interface_equation_a(this, a, g, gt, dg, dgt, r, rt)
      ! This interface is used for the deferred procedure interface_equation_a.
      import material, green, spin, dp

      class(material),          intent(in   ) :: this
      type(green),     pointer, intent(in   ) :: a
      type(spin),               intent(in   ) :: g, gt, dg, dgt
      type(spin),               intent(inout) :: r, rt
    end subroutine
  end interface

  abstract interface
    subroutine interface_equation_b(this, b, g, gt, dg, dgt, r, rt)
      ! This interface is used for the deferred procedure interface_equation_b.
      import material, green, spin, dp

      class(material),          intent(in   ) :: this
      type(green),     pointer, intent(in   ) :: b
      type(spin),               intent(in   ) :: g, gt, dg, dgt
      type(spin),               intent(inout) :: r, rt
    end subroutine
  end interface
contains
  subroutine material_update(this)
    ! This subroutine updates the current estimate for the state of the material by numerically solving the diffusion equation.
    use bvp_m

    class(material), intent(inout) :: this                       ! Material that will be updated
    real(dp)                       :: u(32,size(this%location))  ! Representation of the retarded Green's functions

    complex(dp)                    :: e                          ! Complex energy relative to the Thouless energy
    class(green),          pointer :: a => null()                ! State at this energy at the left  interface
    class(green),          pointer :: b => null()                ! State at this energy at the right interface
    integer                        :: n                          ! Loop variable

    ! Call the prehook method
    call this%update_prehook

    ! Status information
    if (this%information >= 0) then
      write(stdout,'(a)') color_white // ' :: ' // color_none // trim(this%type_string) // '                                     '
    end if

    ! Loop over the discretized energy levels
    do n=1,size(this%energy)
      block
        ! Declare local block variables
        type(bvp_sol) :: sol ! Workspace for the bvp_solver procedures
        integer       :: m   ! Inner loop variable

        ! Status information
        if (this%information >= 0) then
          write(stdout,'(4x,a,1x,i4,1x,a,1x,i4,1x,a,f0.5,a1)',advance='no') &
            '[',n,'/',size(this%energy),']  ϵ = ',this%energy(n), achar(13)
          flush(stdout)
        end if

        ! Convert all states at this energy level to real-valued state vectors
        do m=1,size(this%location)
          u(:,m) = this%greenr(n,m)
        end do

        ! Calculate the complex energy (relative to the Thouless energy)
        e = cmplx(this%energy(n)/this%thouless, this%scattering/this%thouless, kind=dp)

        ! Update the pointers used to evaluate boundary conditions
        if (associated(this%material_a)) then
          a => this%material_a%greenr(n,ubound(this%material_a%greenr,2))
        end if
        if (associated(this%material_b)) then
          b => this%material_b%greenr(n,lbound(this%material_b%greenr,2))
        end if

        ! Initialize bvp_solver
        sol = bvp_init(32, 16, this%location, u, max_num_subintervals=(size(this%location)*this%scaling))

        ! Solve the differential equation
        sol = bvp_solver(sol, ode, bc, method=this%order, error_control=this%control, tol=this%tolerance, trace=this%information)

        ! Use the results to update the state
        call bvp_eval(sol, this%location, u)
        do m=1,size(this%location)
          this%greenr(n,m) = u(:,m)
        end do

        ! Clean up after bvp_solver
        call bvp_terminate(sol)
      end block
    end do

    ! Call the posthook method
    call this%update_posthook
  contains
    subroutine ode(z, u, f)
      ! Definition of the differential equation u'=f(z,u).
      real(dp), intent(in)  :: z
      real(dp), intent(in)  :: u(32)
      real(dp), intent(out) :: f(32)
      type(spin)            :: g, gt, dg, dgt, d2g, d2gt

      ! Extract the Riccati parameters
      g   = u( 1: 8)
      gt  = u( 9:16)
      dg  = u(17:24)
      dgt = u(25:32)

      ! Calculate the second-derivatives of the Riccati parameters
      call this%diffusion_equation(e, z, g, gt, dg, dgt, d2g, d2gt)
       
      ! Pack the results into a state vector
      f( 1: 8) = dg
      f( 9:16) = dgt
      f(17:24) = d2g
      f(25:32) = d2gt
    end subroutine

    subroutine bc(ua, ub, bca, bcb)
      ! Definition of the boundary conditions bca=g(ua) and bcb=g(ub).
      real(dp), intent(in)  :: ua(32)
      real(dp), intent(in)  :: ub(32)
      real(dp), intent(out) :: bca(16)
      real(dp), intent(out) :: bcb(16)

      type(spin)            :: g1, gt1, dg1, dgt1, r1, rt1
      type(spin)            :: g2, gt2, dg2, dgt2, r2, rt2

      ! State at the left end of the material
      g1   = ua( 1: 8)
      gt1  = ua( 9:16)
      dg1  = ua(17:24)
      dgt1 = ua(25:32)

      ! State at the right end of the material
      g2   = ub( 1: 8)
      gt2  = ub( 9:16)
      dg2  = ub(17:24)
      dgt2 = ub(25:32)

      ! Calculate residuals from the boundary conditions
      call this%interface_equation_a(a, g1, gt1, dg1, dgt1, r1, rt1)
      call this%interface_equation_b(b, g2, gt2, dg2, dgt2, r2, rt2)

      ! Pack the results into state vectors
      bca(1: 8) = r1
      bca(9:16) = rt1
      bcb(1: 8) = r2
      bcb(9:16) = rt2
    end subroutine
  end subroutine
end module
