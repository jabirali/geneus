!> Author:   Jabir Ali Ouassou
!> Date:     2015-07-11
!> Category: Foundation
!>
!> This module defines the data type 'propagator', which represents the propagator at a given position and energy. This
!> is done by internally storing the Riccati parameters γ and γ~ and their first derivatives dγ/dz and dγ~/dz. These are
!> sufficient to reconstruct the normal Green's function g and anomalous Green's function f, and derived quantities like
!> the density of states.  To make it easier to interact with differential equation solvers, which often operate on real
!> state vectors, the assignment operator is overloaded so objects can be easily imported/exported to a real vector(32).

module propagator_m
  use math_m
  use spin_m
  implicit none
  private

  ! Public interface
  public propagator, assignment(=)

  ! Type declaration
  type propagator
    type(spin) :: g                                 ! Riccati parameter γ
    type(spin) :: gt                                ! Riccati parameter γ~
    type(spin) :: dg                                ! Derivative dγ /dz
    type(spin) :: dgt                               ! Derivative dγ~/dz
  contains
    procedure  :: matrix    => propagator_matrix    ! Matrix representation of the entire Green's function
    procedure  :: get_g     => propagator_get_g     ! Normal Green's function g
    procedure  :: get_gt    => propagator_get_gt    ! Normal Green's function g~
    procedure  :: get_f     => propagator_get_f     ! Anomal Green's function f
    procedure  :: get_ft    => propagator_get_ft    ! Anomal Green's function f~
    procedure  :: get_f_s   => propagator_get_f_s   ! Singlet component of f
    procedure  :: get_ft_s  => propagator_get_ft_s  ! Singlet component of f~
    procedure  :: get_f_t   => propagator_get_f_t   ! Triplet component of f
    procedure  :: get_ft_t  => propagator_get_ft_t  ! Triplet component of f~
    procedure  :: get_f_ts  => propagator_get_f_ts  ! Short-range triplet component of f
    procedure  :: get_ft_ts => propagator_get_ft_ts ! Short-range triplet component of f~
    procedure  :: get_f_tl  => propagator_get_f_tl  ! Long-range triplet component of f
    procedure  :: get_ft_tl => propagator_get_ft_tl ! Long-range triplet component of f~
    procedure  :: get_dos   => propagator_get_dos   ! Local density of states
    procedure  :: print     => propagator_print     ! Prints the propagator object to standard out
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
    ! as the type(spin) constructors will do it automatically by default
    continue
  end function

  pure function propagator_construct_riccati(g, gt, dg, dgt) result(this)
    !! Construct an arbitrary state by explicitly providing the Riccati parameters.
    type(spin), intent(in)           :: g    !! Riccati parameter
    type(spin), intent(in)           :: gt   !! Riccati parameter (tilde conjugated)
    type(spin), intent(in), optional :: dg   !! Derivative of the Riccati parameter
    type(spin), intent(in), optional :: dgt  !! Derivative of the Riccati parameter (tilde conjugated)
    type(propagator)                 :: this !! Constructed object

    ! Copy Riccati parameters into the new object
    this % g  = g
    this % gt = gt

    ! Copy the derivatives if available
    if (present(dg)) then
      this % dg = dg
    end if
    if (present(dgt)) then
      this % dgt = dgt
    end if
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
  end function

  pure function propagator_matrix(this) result(matrix)
    !! Calculates the 4×4 Green's function matrix from the Riccati parameters of the Green's function object.
    class(propagator), intent(in) :: this        !! Green's function object
    complex(wp)                   :: matrix(4,4) !! Green's function matrix
    type(spin)                    :: N, Nt

    associate(g => this % g, gt => this % gt, I => pauli0, M => matrix)
      ! Calculate the normalization matrices
      N  = spin_inv( I - g*gt )
      Nt = spin_inv( I - gt*g )

      ! Calculate the 4×4 Green's function
      M(1:2,1:2) = (+2.0_wp) * N  - I
      M(1:2,3:4) = (+2.0_wp) * N  * g
      M(3:4,1:2) = (-2.0_wp) * Nt * gt
      M(3:4,3:4) = (-2.0_wp) * Nt + I
    end associate
  end function

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

  pure subroutine propagator_import_rvector(a, b)
    ! Defines assignment from a real vector to a propagator object.
    type(propagator), intent(out) :: a
    real(wp),         intent(in)  :: b(32)

    a%g   = b( 1: 8) 
    a%gt  = b( 9:16) 
    a%dg  = b(17:24) 
    a%dgt = b(25:32) 
  end subroutine

  pure function propagator_get_g(this) result(g)
    ! Calculates the normal Green's function g.
    type(spin)                    :: g
    class(propagator), intent(in) :: this

    g = ( pauli0 - this%g * this%gt ) .divl. ( pauli0 + this%g * this%gt )
  end function

  pure function propagator_get_gt(this) result(gt)
    ! Calculates the tilde-conjugated normal Green's function g~.
    type(spin)                    :: gt
    class(propagator), intent(in) :: this

    gt = ( pauli0 - this%gt * this%g ) .divl. ( pauli0 + this%gt * this%g )
  end function

  pure function propagator_get_f(this) result(f)
    ! Calculates the anomalous Green's function f.
    type(spin)                    :: f
    class(propagator), intent(in) :: this

    f = ( pauli0 - this%g * this%gt ) .divl. ( 2.0_wp * this%g )
  end function

  pure function propagator_get_ft(this) result(ft)
    ! Calculates the tilde-conjugated anomalous Green's function f~.
    type(spin)                    :: ft
    class(propagator), intent(in) :: this

    ft = ( pauli0 - this%gt * this%g ) .divl. ( 2.0_wp * this%gt )
  end function

  pure function propagator_get_f_s(this) result(r)
    ! Calculates the singlet component of the anomalous Green's function f.
    complex(wp)                   :: r
    class(propagator), intent(in) :: this

    type(spin)               :: f
    f = this%get_f()

    r = (f%matrix(1,2) - f%matrix(2,1))/2.0_wp
  end function

  pure function propagator_get_ft_s(this) result(r)
    ! Calculates the singlet component of the tilde-conjugated anomalous Green's function f~.
    complex(wp)                   :: r
    class(propagator), intent(in) :: this

    type(spin)               :: ft
    ft = this%get_ft()

    r = (ft%matrix(1,2) - ft%matrix(2,1))/2.0_wp
  end function

  pure function propagator_get_f_t(this) result(r)
    ! Calculates the triplet component of the anomalous Green's function f.
    complex(wp)                   :: r(3)
    class(propagator), intent(in) :: this

    type(spin)               :: f
    f = this%get_f()

    r = [ (f%matrix(2,2) - f%matrix(1,1))/(2.0_wp,0.0_wp), &
          (f%matrix(1,1) + f%matrix(2,2))/(0.0_wp,2.0_wp), &
          (f%matrix(1,2) + f%matrix(2,1))/(2.0_wp,0.0_wp)  ];
  end function

  pure function propagator_get_ft_t(this) result(r)
    ! Calculates the triplet component of the tilde-conjugated anomalous Green's function f~.
    complex(wp)                   :: r(3)
    class(propagator), intent(in) :: this

    type(spin)               :: ft
    ft = this%get_ft()

    r = [ (ft%matrix(2,2) - ft%matrix(1,1))/(2.0_wp,0.0_wp), &
          (ft%matrix(1,1) + ft%matrix(2,2))/(0.0_wp,2.0_wp), &
          (ft%matrix(1,2) + ft%matrix(2,1))/(2.0_wp,0.0_wp)  ];
  end function

  pure function propagator_get_f_ts(this, h) result(r)
    ! Calculates the short-range triplet component of the anomalous Green's function f,
    ! i.e. the triplet component of f along the magnetic exchange field vector h.
    complex(wp)                   :: r(3)
    real(wp),          intent(in) :: h(3)
    class(propagator), intent(in) :: this

    real(wp)                 :: u(3)
    u = h/(norm2(h)+eps)

    r = dot_product(u,this%get_f_t()) * u
  end function

  pure function propagator_get_ft_ts(this, h) result(r)
    ! Calculates the short-range triplet component of the tilde-conjugated anomalous Green's function f~,
    ! i.e. the triplet component of f~ along the magnetic exchange field h.
    complex(wp)                   :: r(3)
    real(wp),          intent(in) :: h(3)
    class(propagator), intent(in) :: this

    real(wp)                 :: u(3)
    u = h/(norm2(h)+eps)

    r = dot_product(u,this%get_ft_t()) * u
  end function

  pure function propagator_get_f_tl(this, h) result(r)
    ! Calculates the long-range triplet component of the anomalous Green's function f,
    ! i.e. the triplet component of f perpendicular to the magnetic exchange field h.
    complex(wp)                   :: r(3)
    real(wp),          intent(in) :: h(3)
    class(propagator), intent(in) :: this

    r = this%get_f_t() - this%get_f_ts(h)
  end function

  pure function propagator_get_ft_tl(this, h) result(r)
    ! Calculates the long-range triplet component of the tilde-conjugated anomalous Green's function f~,
    ! i.e. the triplet component of f~ perpendicular to the magnetic exchange field h.
    complex(wp)                   :: r(3)
    real(wp),          intent(in) :: h(3)
    class(propagator), intent(in) :: this

    r = this%get_ft_t() - this%get_ft_ts(h)
  end function

  pure function propagator_get_dos(this) result(ldos)
    ! Calculates the local density of states.
    real(wp)                      :: ldos
    class(propagator), intent(in) :: this

    type(spin)               :: g
    g = this%get_g()

    ldos = 0.5_wp*re(g%trace())
  end function

  impure subroutine propagator_print(this)
    ! Prints the propagator object to stdout.
    class(propagator), intent(in) :: this 

    ! Print the matrix elements
    call this%g%print('Riccati parameter gamma ')
    call this%gt%print('Riccati parameter gamma~')
    call this%dg%print('Derivative dgamma / dz')
    call this%dgt%print('Derivative dgamma~/ dz')
  end subroutine
end module
