!> Author:   Jabir Ali Ouassou
!> Date:     2015-07-11
!> Category: Foundation
!>
!> This module defines the data type 'propagator', which represents the propagator at a given position and energy. This
!> is done by internally storing the Riccati parameters γ and γ~ and their first derivatives dγ/dz and dγ~/dz. These are
!> sufficient to reconstruct the 4×4 matrix representation of the propagator, and associated physical quantities such as
!> the density of states.  To make it easier to interact with differential equation solvers, which often operate on real
!> state vectors, the assignment operator is overloaded so that objects can be transformed to 32-element real arrays.

module propagator_m
  use :: math_m
  use :: spin_m
  use :: nambu_m
  private

  ! Public interface
  public propagator, assignment(=)

  ! Type declaration
  type propagator
    type(spin) :: g                                 !! Riccati parameter γ
    type(spin) :: gt                                !! Riccati parameter γ~
    type(spin) :: dg                                !! Derivative dγ /dz
    type(spin) :: dgt                               !! Derivative dγ~/dz
    type(spin) :: N                                 !! Normalization N  = (1 - γγ~)^-1
    type(spin) :: Nt                                !! Normalization Nt = (1 - γ~γ)^-1
  contains
    ! Accessors for the propagator matrices represented by this object
    procedure  :: retarded           => propagator_retarded           !! Retarded propagator (G^R)
    procedure  :: advanced           => propagator_advanced           !! Advanced propagator (G^A)
    procedure  :: retarded_gradient  => propagator_retarded_gradient  !! Retarded propagator gradient (dG^R/dz)
    procedure  :: advanced_gradient  => propagator_advanced_gradient  !! Advanced propagator gradient (dG^A/dz)

    ! Accessors for physical quantities that derive from this object
    procedure  :: decompose          => propagator_decompose          !! Singlet/triplet decomposition of the anomalous retarded propagator (f)
    procedure  :: density            => propagator_density            !! Density of states at this position z and energy ε
    procedure  :: current            => propagator_current            !! Spectral current  at this position z and energy ε
  end type

  ! Type constructor
  interface propagator
    module procedure propagator_construct_zero, propagator_construct_riccati, propagator_construct_bcs
  end interface

  ! Assignment operator
  interface assignment(=)
    module procedure propagator_import_rvector, propagator_export_rvector, propagator_export_cmatrix
  end interface
contains
  pure function propagator_construct_zero() result(this)
    ! Constructs a state corresponding to a normal metal, which has all the Riccati parameters set to zero.
    type(propagator) :: this

    ! There is no need to explicitly set the Riccati parameters to zero, 
    ! as the type(spin) constructors will do it automatically by default.

    ! Update normalization matrices
    this % N  = pauli0
    this % Nt = pauli0
  end function

  pure function propagator_construct_riccati(g, gt) result(this)
    !! Construct an arbitrary state by explicitly providing the Riccati parameters.
    type(spin), intent(in)           :: g    !! Riccati parameter
    type(spin), intent(in)           :: gt   !! Riccati parameter (tilde conjugated)
    type(propagator)                 :: this !! Constructed object

    ! Copy Riccati parameters into the new object
    this % g  = g
    this % gt = gt

    ! Update normalization matrices
    this % N  = inverse( pauli0 - g*gt )
    this % Nt = inverse( pauli0 - gt*g )
  end function

  pure function propagator_construct_bcs(energy, gap) result(this)
    ! Constructs a state corresponding to a BCS superconductor at some given energy, which may have an imaginary
    ! term representing inelastic scattering. The second argument 'gap' is the superconducting order parameter Δ.
    type(propagator)        :: this      ! Propagator object that will be constructed
    complex(wp), intent(in) :: energy    ! Quasiparticle energy (including inelastic scattering contribution)
    complex(wp), intent(in) :: gap       ! Superconducting order parameter (including superconducting phase)

    real(wp)                :: p
    complex(wp)             :: t, u
    complex(wp)             :: a, b

    ! Calculate the superconducting gap and phase
    u = abs(gap)/energy
    p = atan2(im(gap), re(gap))

    ! Calculate the θ-parameter
    t = (log(1+u)-log(1-u))/2

    ! Calculate the scalar Riccati parameters a and b
    a =  (exp(+t)-exp(-t))/(2+exp(+t)+exp(-t)) * exp( (0,+1) * p )
    b = -(exp(+t)-exp(-t))/(2+exp(+t)+exp(-t)) * exp( (0,-1) * p )

    ! Calculate the matrix Riccati parameters γ and γ~
    this % g  = a * ((0.0_wp,1.0_wp) * pauli2)
    this % gt = b * ((0.0_wp,1.0_wp) * pauli2)

    ! Update normalization matrices
    this % N  = inverse( pauli0 - this%g  * this%gt )
    this % Nt = inverse( pauli0 - this%gt * this%g  )
  end function

  pure function propagator_retarded(this) result(r)
    !! Calculates the 4×4 retarded propagator (G^R).
    class(propagator), intent(in) :: this   !! Propagator object
    type(nambu)                   :: r      !! Propagator matrix

    associate(g => this % g, gt => this % gt, &
              N => this % N, Nt => this % Nt, &
              I => pauli0,   M  => r % matrix )
      M(1:2,1:2) = (+2.0_wp) * N  - I
      M(1:2,3:4) = (+2.0_wp) * N  * g
      M(3:4,1:2) = (-2.0_wp) * Nt * gt
      M(3:4,3:4) = (-2.0_wp) * Nt + I
    end associate
  end function

  pure function propagator_retarded_gradient(this) result(r)
    !! Calculates the 4×4 retarded propagator gradient (dG^R/dz).
    class(propagator), intent(in) :: this   !! Propagator object
    type(nambu)                   :: r      !! Propagator gradient

    associate(g  => this % g,  gt  => this % gt,  &
              dg => this % dg, dgt => this % dgt, &
              N  => this % N,  Nt  => this % Nt,  &
              I  => pauli0,    M   => r % matrix  )
      M(1:2,1:2) = (+1.0_wp) * N  * (dg*gt  +  g*dgt) * N
      M(1:2,3:4) = (+2.0_wp) * N  * (dg  - g *dgt*g ) * Nt
      M(3:4,1:2) = (-2.0_wp) * Nt * (dgt - gt*dg *gt) * N
      M(3:4,3:4) = (-1.0_wp) * Nt * (dgt*g  +  gt*dg) * Nt
    end associate
  end function

  pure function propagator_advanced(this) result(r)
    !! Calculates the 4×4 advanced propagator (G^A).
    class(propagator), intent(in) :: this   !! Propagator object
    type(nambu)                   :: r      !! Propagator matrix

    r = this % retarded()
    r = nambuv(4) * transpose(conjg(-r % matrix)) * nambuv(4)
  end function

  pure function propagator_advanced_gradient(this) result(r)
    !! Calculates the 4×4 advanced propagator gradient (dG^A/dz).
    class(propagator), intent(in) :: this   !! Propagator object
    type(nambu)                   :: r      !! Propagator gradient

    r = this % retarded_gradient()
    r = nambuv(4) * transpose(conjg(-r % matrix)) * nambuv(4)
  end function

  pure subroutine propagator_import_rvector(a, b)
    ! Defines assignment from a real vector to a propagator object.
    type(propagator), intent(out) :: a
    real(wp),         intent(in)  :: b(32)

    a%g   = b( 1: 8) 
    a%gt  = b( 9:16) 
    a%dg  = b(17:24) 
    a%dgt = b(25:32) 

    ! Update normalization matrices
    a % N  = inverse( pauli0 - a%g  * a%gt )
    a % Nt = inverse( pauli0 - a%gt * a%g  )
  end subroutine

  pure subroutine propagator_export_cmatrix(a, b)
    ! Defines assignment from a propagator object to a complex matrix.
    complex(wp),      intent(out) :: a(4,4)
    type(propagator), intent(in)  :: b

    a = b % retarded()
  end subroutine

  pure subroutine propagator_export_rvector(a, b)
    ! Defines assignment from a propagator object to a real vector.
    real(wp),         intent(out) :: a(32)
    type(propagator), intent(in)  :: b

    a( 1: 8) = b%g
    a( 9:16) = b%gt
    a(17:24) = b%dg
    a(25:32) = b%dgt
  end subroutine

  pure subroutine propagator_decompose(this, f, ft, df, dft)
    !! Performs a singlet/triplet decomposition of the anomalous retarded propagators (f, f~)
    !! and their derivatives (df/dz, df~/dz). Each of these are returned as a complex 4-vector,
    !! where f(0) corresponds to the singlet component and f(1:3) to the triplet component.
    class(propagator),                     intent(in)  :: this  !! Propagator object
    complex(wp), dimension(0:3), optional, intent(out) :: f     !! Anomalous propagator (f )
    complex(wp), dimension(0:3), optional, intent(out) :: ft    !! Anomalous propagator (f~)
    complex(wp), dimension(0:3), optional, intent(out) :: df    !! Anomalous propagator gradient (df /dz)
    complex(wp), dimension(0:3), optional, intent(out) :: dft   !! Anomalous propagator gradient (df~/dz)
    type(nambu)                                        :: G     !  4×4 retarded propagator (G^R)

    ! Calculate and decompose the anomalous propagator
    if (present(f) .or. present(ft)) then
      ! Extract the retarded propagator
      G  = this % retarded()

      ! Perform the singlet/triplet decomposition
      if (present(f )) then
        f   = trace(pauli * ( G % matrix(1:2, 3:4) * pauli2)) * (0.0_wp,-0.5_wp)
      end if
      if (present(ft)) then
        ft  = trace(pauli * ( G % matrix(3:4, 1:2) * pauli2)) * (0.0_wp,+0.5_wp)
      end if
    end if

    ! Calculate and decompose the gradient
    if (present(df) .or. present(dft)) then
      ! Extract the retarded propagator gradient
      G  = this % retarded_gradient()

      ! Perform the singlet/triplet decomposition
      if (present(df )) then
        df   = trace(pauli * ( G % matrix(1:2, 3:4) * pauli2)) * (0.0_wp,-0.5_wp)
      end if
      if (present(dft)) then
        dft  = trace(pauli * ( G % matrix(3:4, 1:2) * pauli2)) * (0.0_wp,+0.5_wp)
      end if
    end if
  end subroutine

  pure function propagator_density(this) result(r)
    ! Calculates the local density of states.
    class(propagator), intent(in) :: this
    real(wp)                      :: r

    r = re(trace(this % N)) - 1
  end function

  pure function propagator_current(this) result(r)
    !! Calculates the spectral current at zero temperature. The result is a 4-vector,
    !! where element 0 is the charge current, and elements 1:3 are the spin currents.
    !!
    !! @TODO: 
    !!   The equation below should be changed from:
    !!     r = 8 * re(trace(pauli * k))
    !!   to the more logical form:
    !!     r = 2 * re(trace(pauli * k))
    !!   This would change the current unit to J₀=eN₀∆₀²ξ²A/L, without a factor 1/4.
    class(propagator), intent(in) :: this
    real(wp)                      :: r(0:3)
    type(spin)                    :: k

    associate(g  => this % g,  dg  => this % dg,  N  => this % N, &
              gt => this % gt, dgt => this % dgt, Nt => this % Nt )
      k = N*(dg*gt-g*dgt)*N - conjg(Nt*(dgt*g-gt*dg)*Nt)
      r = 8 * re(trace(pauli * k))
    end associate
  end function
end module
