! This module defines the data type 'green', which represents the Green's function at a given position and energy. This
! is done by internally storing the Riccati parameters γ and γ~ and their first derivatives dγ/dx and dγ~/dx. These are
! sufficient to reconstruct the normal Green's function g and anomalous Green's function f, and derived quantities like
! the density of states.  To make it easier to interact with differential equation solvers, which often operate on real
! state vectors, the assignment operator is overloaded so 'green' objects can be imported from and exported to real(32).
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-11
! Updated: 2015-07-20

module mod_green
  use mod_system
  use mod_spin
  implicit none

  ! Type declaration
  type green
    type(spin) :: g                            ! Riccati parameter γ  
    type(spin) :: gt                           ! Riccati parameter γ~
    type(spin) :: dg                           ! Derivative dγ /dz
    type(spin) :: dgt                          ! Derivative dγ~/dz
    contains
    procedure  :: get_g     => green_get_g     ! Normal Green's function g
    procedure  :: get_gt    => green_get_gt    ! Normal Green's function g~
    procedure  :: get_f     => green_get_f     ! Anomal Green's function f
    procedure  :: get_ft    => green_get_ft    ! Anomal Green's function f~
    procedure  :: get_f_s   => green_get_f_s   ! Singlet component of f
    procedure  :: get_ft_s  => green_get_ft_s  ! Singlet component of f~
    procedure  :: get_f_t   => green_get_f_t   ! Triplet component of f
    procedure  :: get_ft_t  => green_get_ft_t  ! Triplet component of f~
    procedure  :: get_f_ts  => green_get_f_ts  ! Short-range triplet component of f
    procedure  :: get_ft_ts => green_get_ft_ts ! Short-range triplet component of f~
    procedure  :: get_f_tl  => green_get_f_tl  ! Long-range triplet component of f
    procedure  :: get_ft_tl => green_get_ft_tl ! Long-range triplet component of f~
    procedure  :: get_dos   => green_get_dos   ! Local density of states
    procedure  :: print     => green_print     ! Prints the green object to standard out
  end type

  ! Type constructor
  interface green
    module procedure green_construct_bcs
  end interface

  ! Assignment operator
  interface assignment(=)
    module procedure green_import_rvector, green_export_rvector
  end interface
contains
  pure function green_construct_bcs(energy, gap) result(this)
    ! Constructs a state corresponding to a BCS superconductor at some given energy, which may have an imaginary
    ! term representing inelastic scattering. The second argument 'gap' is the superconducting order parameter Δ.
    type(green)             :: this      ! Green's function object that will be constructed
    complex(dp), intent(in) :: energy    ! Quasiparticle energy (including inelastic scattering contribution)
    complex(dp), intent(in) :: gap       ! Superconducting order parameter (including superconducting phase)

    ! Calculate the θ-parameter (t), superconducting phase (p), and scalar Riccati parameters (a,b)
    complex(dp)             :: t
    real(dp)                :: p
    complex(dp)             :: a, b

    t = atanh(abs(gap)/energy)
    p = atan(aimag(gap)/dble(gap))
    a =  sinh(t)/(1+cosh(t)) * exp( (0, 1) * p )
    b = -sinh(t)/(1+cosh(t)) * exp( (0,-1) * p )

    ! Calculate the matrix Riccati parameters γ and γ~
    this%g  = [(0.0_dp,0.0_dp), a, -a, (0.0_dp,0.0_dp)]
    this%gt = [(0.0_dp,0.0_dp), b, -b, (0.0_dp,0.0_dp)]
  end function

  pure subroutine green_export_rvector(a, b)
    ! Defines assignment from a green object to a real vector.
    real(dp),    intent(out) :: a(32)
    type(green), intent(in)  :: b

    a( 1: 8) = b%g
    a( 9:16) = b%gt
    a(17:24) = b%dg
    a(25:32) = b%dgt
  end subroutine

  pure subroutine green_import_rvector(a, b)
    ! Defines assignment from a real vector to a green object.
    type(green), intent(out) :: a
    real(dp),    intent(in)  :: b(32)

    a%g   = b( 1: 8) 
    a%gt  = b( 9:16) 
    a%dg  = b(17:24) 
    a%dgt = b(25:32) 
  end subroutine

  function green_get_g(this) result(g)
    ! Calculates the normal Green's function g.
    type(spin)               :: g
    class(green), intent(in) :: this

    g = ( pauli0 - this%g * this%gt ) .divl. ( pauli0 + this%g * this%gt )
  end function

  function green_get_gt(this) result(gt)
    ! Calculates the tilde-conjugated normal Green's function g~.
    type(spin)               :: gt
    class(green), intent(in) :: this

    gt = ( pauli0 - this%gt * this%g ) .divl. ( pauli0 + this%gt * this%g )
  end function

  function green_get_f(this) result(f)
    ! Calculates the anomalous Green's function f.
    type(spin)               :: f
    class(green), intent(in) :: this

    f = ( pauli0 - this%g * this%gt ) .divl. ( 2.0_dp * this%g )
  end function

  function green_get_ft(this) result(ft)
    ! Calculates the tilde-conjugated anomalous Green's function f~.
    type(spin)               :: ft
    class(green), intent(in) :: this

    ft = ( pauli0 - this%gt * this%g ) .divl. ( 2.0_dp * this%gt )
  end function

  function green_get_f_s(this) result(r)
    ! Calculates the singlet component of the anomalous Green's function f.
    complex(dp)              :: r
    class(green), intent(in) :: this
    
    type(spin)               :: f
    f = this%get_f()

    r = (f%matrix(1,2) - f%matrix(2,1))/2.0_dp
  end function

  function green_get_ft_s(this) result(r)
    ! Calculates the singlet component of the tilde-conjugated anomalous Green's function f~.
    complex(dp)              :: r
    class(green), intent(in) :: this
    
    type(spin)               :: ft
    ft = this%get_ft()

    r = (ft%matrix(1,2) - ft%matrix(2,1))/2.0_dp
  end function

  function green_get_f_t(this) result(r)
    ! Calculates the triplet component of the anomalous Green's function f.
    complex(dp)              :: r(3)
    class(green), intent(in) :: this

    type(spin)               :: f
    f = this%get_f()

    r = [ (f%matrix(2,2) - f%matrix(1,1))/(2.0_dp,0.0_dp), &
          (f%matrix(1,1) + f%matrix(2,2))/(0.0_dp,2.0_dp), &
          (f%matrix(1,2) + f%matrix(2,1))/(2.0_dp,0.0_dp)  ];
  end function

  function green_get_ft_t(this) result(r)
    ! Calculates the triplet component of the tilde-conjugated anomalous Green's function f~.
    complex(dp)              :: r(3)
    class(green), intent(in) :: this

    type(spin)               :: ft
    ft = this%get_ft()

    r = [ (ft%matrix(2,2) - ft%matrix(1,1))/(2.0_dp,0.0_dp), &
          (ft%matrix(1,1) + ft%matrix(2,2))/(0.0_dp,2.0_dp), &
          (ft%matrix(1,2) + ft%matrix(2,1))/(2.0_dp,0.0_dp)  ];
  end function

  function green_get_f_ts(this, h) result(r)
    ! Calculates the short-range triplet component of the anomalous Green's function f,
    ! i.e. the triplet component of f along the magnetic exchange field vector h.
    complex(dp)              :: r(3)
    real(dp),     intent(in) :: h(3)
    class(green), intent(in) :: this

    real(dp)                 :: u(3)
    u = h/(norm2(h)+1e-16_dp)

    r = dot_product(u,this%get_f_t()) * u
  end function

  function green_get_ft_ts(this, h) result(r)
    ! Calculates the short-range triplet component of the tilde-conjugated anomalous Green's function f~,
    ! i.e. the triplet component of f~ along the magnetic exchange field h.
    complex(dp)              :: r(3)
    real(dp),     intent(in) :: h(3)
    class(green), intent(in) :: this

    real(dp)                 :: u(3)
    u = h/(norm2(h)+1e-16_dp)

    r = dot_product(u,this%get_ft_t()) * u
  end function

  function green_get_f_tl(this, h) result(r)
    ! Calculates the long-range triplet component of the anomalous Green's function f,
    ! i.e. the triplet component of f perpendicular to the magnetic exchange field h.
    complex(dp)              :: r(3)
    real(dp),     intent(in) :: h(3)
    class(green), intent(in) :: this

    r = this%get_f_t() - this%get_f_ts(h)
  end function

  function green_get_ft_tl(this, h) result(r)
    ! Calculates the long-range triplet component of the tilde-conjugated anomalous Green's function f~,
    ! i.e. the triplet component of f~ perpendicular to the magnetic exchange field h.
    complex(dp)              :: r(3)
    real(dp),     intent(in) :: h(3)
    class(green), intent(in) :: this

    r = this%get_ft_t() - this%get_ft_ts(h)
  end function

  function green_get_dos(this) result(ldos)
    ! Calculates the local density of states.
    real(dp)                 :: ldos
    class(green), intent(in) :: this

    type(spin)               :: g
    g = this%get_g()

    ldos = real(g%trace(), kind=dp)/2.0_dp
  end function

  subroutine green_print(this)
    ! Prints the green object to stdout.
    class(green), intent(in) :: this 

    ! Print the matrix elements
    call this%g%print('Riccati parameter γ ')
    call this%gt%print('Riccati parameter γ~')
    call this%dg%print('Derivative dγ / dz')
    call this%dgt%print('Derivative dγ~/ dz')
  end subroutine
end module
