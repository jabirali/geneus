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
  use math_m
  use spin_m
  implicit none
  private

  ! Public interface
  public propagator, assignment(=)

  ! Type declaration
  type propagator
    type(spin) :: g                                 !! Riccati parameter γ
    type(spin) :: gt                                !! Riccati parameter γ~
    type(spin) :: dg                                !! Derivative dγ /dz
    type(spin) :: dgt                               !! Derivative dγ~/dz
    type(spin) :: N                                 !! Normalization N  = inv(I - γγ~)
    type(spin) :: Nt                                !! Normalization Nt = inv(I - γ~γ)
  contains
    ! Accessors for propagator components
    procedure  :: matrix    => propagator_matrix    !! 4×4 matrix representation of the propagator
    procedure  :: singlet   => propagator_singlet   !! Singlet component of the anomalous propagator
    procedure  :: singlett  => propagator_singlett  !! Singlet component of the anomalous propagator (tilde-conjugated)
    procedure  :: triplet   => propagator_triplet   !! Triplet component of the anomalous propagator
    procedure  :: triplett  => propagator_triplett  !! Triplet component of the anomalous propagator (tilde-conjugated)

    ! Accessors for derived physical quantities
    procedure  :: density   => propagator_density   !! Density of states at this position and energy
    procedure  :: current   => propagator_current   !! Spectral current  at this position and energy

    ! Miscellaneous procedures
    procedure  :: print     => propagator_print     !! Prints the propagator object to standard out
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
    this % N  = spin_inv( pauli0 - g*gt )
    this % Nt = spin_inv( pauli0 - gt*g )
  end function

  pure function propagator_construct_bcs(energy, gap) result(this)
    ! Constructs a state corresponding to a BCS superconductor at some given energy, which may have an imaginary
    ! term representing inelastic scattering. The second argument 'gap' is the superconducting order parameter Δ.
    type(propagator)        :: this      ! Green's function object that will be constructed
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
    this % g  = [(0.0_wp,0.0_wp), a, -a, (0.0_wp,0.0_wp)]
    this % gt = [(0.0_wp,0.0_wp), b, -b, (0.0_wp,0.0_wp)]

    ! Update normalization matrices
    this % N  = spin_inv( pauli0 - this%g  * this%gt )
    this % Nt = spin_inv( pauli0 - this%gt * this%g  )
  end function

  pure function propagator_matrix(this) result(matrix)
    !! Calculates the 4×4 Green's function matrix from the Riccati parameters of the Green's function object.
    class(propagator), intent(in) :: this        !! Green's function object
    complex(wp)                   :: matrix(4,4) !! Green's function matrix

    associate(g => this % g, gt => this % gt, N => this % N, Nt => this % Nt, I => pauli0, M => matrix)
      M(1:2,1:2) = (+2.0_wp) * N  - I
      M(1:2,3:4) = (+2.0_wp) * N  * g
      M(3:4,1:2) = (-2.0_wp) * Nt * gt
      M(3:4,3:4) = (-2.0_wp) * Nt + I
    end associate
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
    a % N  = spin_inv( pauli0 - a%g  * a%gt )
    a % Nt = spin_inv( pauli0 - a%gt * a%g  )
  end subroutine

  pure subroutine propagator_export_cmatrix(a, b)
    ! Defines assignment from a propagator object to a complex matrix.
    complex(wp),      intent(out) :: a(4,4)
    type(propagator), intent(in)  :: b

    a = b % matrix()
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

  pure function propagator_singlet(this) result(r)
    ! Calculates the singlet component of the anomalous Green's function f.
    class(propagator), intent(in) :: this
    complex(wp)                   :: r
    type(spin)                    :: f

    f = this % N * this % g
    r = f % matrix(1,2) - f % matrix(2,1)
  end function

  pure function propagator_singlett(this) result(r)
    ! Calculates the singlet component of the tilde-conjugated anomalous Green's function f~.
    class(propagator), intent(in) :: this
    complex(wp)                   :: r
    type(spin)                    :: ft

    ft = this % Nt * this % gt
    r  = ft % matrix(1,2) - ft % matrix(2,1)
  end function

  pure function propagator_triplet(this) result(r)
    ! Calculates the triplet component of the anomalous Green's function f.
    class(propagator), intent(in) :: this
    complex(wp)                   :: r(3)
    type(spin)                    :: f

    f = this % N * this % g

    r = [ (f%matrix(2,2) - f%matrix(1,1))*(1.0_wp, 0.0_wp), &
          (f%matrix(1,1) + f%matrix(2,2))*(0.0_wp,-1.0_wp), &
          (f%matrix(1,2) + f%matrix(2,1))*(1.0_wp, 0.0_wp)  ];
  end function

  pure function propagator_triplett(this) result(r)
    ! Calculates the triplet component of the tilde-conjugated anomalous Green's function f~.
    class(propagator), intent(in) :: this
    complex(wp)                   :: r(3)
    type(spin)                    :: ft

    ft = this % Nt * this % gt

    r = [ (ft%matrix(2,2) - ft%matrix(1,1))*(1.0_wp, 0.0_wp), &
          (ft%matrix(1,1) + ft%matrix(2,2))*(0.0_wp,-1.0_wp), &
          (ft%matrix(1,2) + ft%matrix(2,1))*(1.0_wp, 0.0_wp)  ];
  end function


  pure function propagator_density(this) result(r)
    ! Calculates the local density of states.
    class(propagator), intent(in) :: this
    real(wp)                      :: r

    r = re(spin_trace(this % N)) - 1
  end function

  pure function propagator_current(this) result(r)
    !! Calculates the spectral current at zero temperature. The result is a 4-vector,
    !! where element 0 is the charge current, and elements 1:3 are the spin currents.
    class(propagator), intent(in) :: this
    real(wp)                      :: r(0:3)
    type(spin)                    :: p(0:3)
    type(spin)                    :: k

    associate(g  => this % g,  dg  => this % dg,  N =>  this % N, &
              gt => this % gt, dgt => this % dgt, Nt => this % Nt )
      p = 8.0_wp * [pauli0, pauli1, pauli2, pauli3]
      k = N*(dg*gt-g*dgt)*N - conjg(Nt*(dgt*g-gt*dg)*Nt)
      r = re(spin_trace(p*k))
    end associate
  end function

  impure subroutine propagator_print(this)
    ! Prints the propagator object to stdout.
    class(propagator), intent(in) :: this 

    ! Print the matrix elements
    call this % g   % print('Riccati parameter γ ')
    call this % gt  % print('Riccati parameter γ~')
    call this % dg  % print('Derivative dγ / dz')
    call this % dgt % print('Derivative dγ~/ dz')
  end subroutine
end module
