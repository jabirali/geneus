!> Author:   Jabir Ali Ouassou
!> Category: Foundation
!>
!> This module defines the data type 'propagator', which represents the propagator at a given position and energy. 
!> The equilibrium propagators (retarded and advanced) are stored internally using the Riccati parameters γ and γ~
!> and their derivatives, while the nonequilibrium propagator (Keldysh) is represented by taking the traces of the
!> distribution function and its derivatives. These quantities are together sufficient to reconstruct the full 8×8
!> propagator and its derivatives, and can be used to calculate associated physical quantities such as the density
!> of states, charge currents, spin currents, heat currents, spin-heat currents, and various accumulation effects.

module propagator_m
  use :: math_m
  use :: spin_m
  use :: nambu_m
  private

  ! Public interface
  public propagator, assignment(=)

  ! Type declaration
  type propagator
    ! Riccati parametrization of equilibrium propagators (retarded and advanced)
    type(spin) :: g                                                             !! Riccati parameter γ
    type(spin) :: gt                                                            !! Riccati parameter γ~
    type(spin) :: dg                                                            !! Riccati parameter ∇γ
    type(spin) :: dgt                                                           !! Riccati parameter ∇γ~
    type(spin) :: d2g                                                           !! Riccati parameter ∇²γ
    type(spin) :: d2gt                                                          !! Riccati parameter ∇²γ~
    type(spin) :: N                                                             !! Riccati normalization N
    type(spin) :: Nt                                                            !! Riccati normalization N~

    ! Distribution-trace parametrization of nonequilibrium propagators (Keldysh)
    real(wp), dimension(0:7) :: h   = [1,0,0,0,0,0,0,0]                         !! Distribution-trace H
    real(wp), dimension(0:7) :: dh  = [0,0,0,0,0,0,0,0]                         !! Distribution-trace ∇H
    real(wp), dimension(0:7) :: d2h = [0,0,0,0,0,0,0,0]                         !! Distribution-trace ∇²H
  contains
    ! Accessors for the propagator matrices represented by this object
    procedure  :: retarded                => propagator_retarded                !! Retarded propagator G^R
    procedure  :: retarded_gradient       => propagator_retarded_gradient       !! Retarded propagator ∇G^R
  ! procedure  :: retarded_laplacian      => propagator_retarded_laplacian      !! Retarded propagator ∇²G^R

    procedure  :: advanced                => propagator_advanced                !! Advanced propagator G^A
    procedure  :: advanced_gradient       => propagator_advanced_gradient       !! Advanced propagator ∇G^A
  ! procedure  :: advanced_laplacian      => propagator_advanced_laplacian      !! Advanced propagator ∇²G^A

    procedure  :: keldysh                 => propagator_keldysh                 !! Keldysh propagator G^K
    procedure  :: keldysh_gradient        => propagator_keldysh_gradient        !! Keldysh propagator ∇G^K
  ! procedure  :: keldysh_laplacian       => propagator_keldysh_laplacian       !! Keldysh propagator ∇²G^K

    procedure  :: distribution            => propagator_distribution            !! Distribution matrix H
    procedure  :: distribution_gradient   => propagator_distribution_gradient   !! Distribution matrix ∇H
  ! procedure  :: distribution_laplacian  => propagator_distribution_laplacian  !! Distribution matrix ∇²H

    ! Accessors for physical quantities that derive from the propagators
    procedure  :: supercurrent       => propagator_supercurrent                 !! Spectral supercurrents
    procedure  :: lossycurrent       => propagator_lossycurrent                 !! Spectral dissipative currents
    procedure  :: accumulation       => propagator_accumulation                 !! Spectral accumulations
    procedure  :: correlation        => propagator_correlation                  !! Spectral correlations
    procedure  :: density            => propagator_density                      !! Spin-resolved density of states
    procedure  :: decompose          => propagator_decompose                    !! Singlet/triplet decomposition
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
    !! Constructs a state corresponding to a normal metal, which has all 
    !! the Riccati parameters set to zero without proximity effects. The
    !! distribution function defaults to equilibrium at zero temperature.
    type(propagator) :: this                !! Constructed object

    ! There is no need to explicitly set the Riccati parameters to zero, 
    ! as the type(spin) constructors will do it automatically by default.

    ! Update normalization matrices
    this % N  = pauli0
    this % Nt = pauli0
  end function

  pure function propagator_construct_riccati(g, gt) result(this)
    !! Construct an arbitrary state by explicitly providing the Riccati parameters.
    !! The distribution function defaults to equilibrium at zero temperature.
    type(propagator)                 :: this !! Constructed object
    type(spin), intent(in)           :: g    !! Riccati parameter γ
    type(spin), intent(in)           :: gt   !! Riccati parameter γ~

    ! Copy Riccati parameters into the new object
    this % g  = g
    this % gt = gt

    ! Update the normalization matrices
    this % N  = inverse( pauli0 - g*gt )
    this % Nt = inverse( pauli0 - gt*g )
  end function

  pure function propagator_construct_bcs(energy, gap) result(this)
    !! Constructs a state corresponding to a BCS superconductor at some given energy,
    !! which may have an imaginary term representing inelastic scattering. The second
    !! argument 'gap' is used to provide the superconducting order parameter Δ.
    !! The distribution function defaults to equilibrium at zero temperature.
    type(propagator)        :: this      !! Constructed object
    complex(wp), intent(in) :: energy    !! Quasiparticle energy (including inelastic scattering contribution)
    complex(wp), intent(in) :: gap       !! Superconducting order parameter (including superconducting phase)

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

    ! Update the normalization matrices
    this % N  = inverse( pauli0 - this%g  * this%gt )
    this % Nt = inverse( pauli0 - this%gt * this%g  )
  end function

  pure function propagator_retarded(this) result(GR)
    !! Calculates the 4×4 retarded propagator G^R.
    class(propagator), intent(in) :: this   !! Propagator object
    type(nambu)                   :: GR     !! Retarded propagator

    ! Construct the propagator from the Riccati parameters
    associate(g => this % g, gt => this % gt, &
              N => this % N, Nt => this % Nt, &
              I => pauli0,   M  => GR % matrix)
      M(1:2,1:2) = (+2.0_wp) * N  - I
      M(1:2,3:4) = (+2.0_wp) * N  * g
      M(3:4,1:2) = (-2.0_wp) * Nt * gt
      M(3:4,3:4) = (-2.0_wp) * Nt + I
    end associate
  end function

  pure function propagator_retarded_gradient(this, gauge) result(dGR)
    !! Calculates the 4×4 retarded propagator gradient ∇G^R. If an optional
    !! gauge field is specified, it returns the gauge-covariant gradient.
    class(propagator),     intent(in) :: this   !! Propagator object
    type(nambu), optional, intent(in) :: gauge  !! Optional gauge field
    type(nambu)                       :: dGR    !! Retarded propagator gradient

    ! Construct the propagator from the Riccati parameters
    associate(g  => this % g,  gt  => this % gt,  &
              dg => this % dg, dgt => this % dgt, &
              N  => this % N,  Nt  => this % Nt,  &
              I  => pauli0,    M   => dGR % matrix)
      M(1:2,1:2) = (+2.0_wp) * N  * (dg*gt  +  g*dgt) * N
      M(1:2,3:4) = (+2.0_wp) * N  * (dg  + g *dgt*g ) * Nt
      M(3:4,1:2) = (-2.0_wp) * Nt * (dgt + gt*dg *gt) * N
      M(3:4,3:4) = (-2.0_wp) * Nt * (dgt*g  +  gt*dg) * Nt
    end associate

    ! Construct the gauge-covariant terms
    if (present(gauge)) then
      associate (A  => gauge, GR => this % retarded(), i  => (0.0_wp,1.0_wp))
        dGR = dGR - i*(A*GR - GR*A)
      end associate
    end if
  end function

  pure function propagator_advanced(this) result(GA)
    !! Calculates the 4×4 advanced propagator G^A.
    class(propagator), intent(in) :: this   !! Propagator object
    type(nambu)                   :: GA     !! Advanced propagator
    type(nambu)                   :: GR     !! Retarded propagator

    ! Calculate the retarded propagator
    GR = this % retarded()

    ! Use the identity GA = -τ₃GR†τ₃
    GA = nambuv(4) * transpose(conjg(-GR % matrix)) * nambuv(4)
  end function

  pure function propagator_advanced_gradient(this, gauge) result(dGA)
    !! Calculates the 4×4 advanced propagator gradient ∇G^A. If an optional
    !! gauge field is specified, it returns the gauge-covariant gradient.
    class(propagator),     intent(in) :: this   !! Propagator object
    type(nambu), optional, intent(in) :: gauge  !! Optional gauge field
    type(nambu)                       :: dGA    !! Advanced propagator gradient
    type(nambu)                       :: dGR    !! Retarded propagator gradient

    ! Calculate the retarded propagator gradient
    dGR = this % retarded_gradient(gauge)

    ! Use the identity GA = -τ₃GR†τ₃
    dGA = nambuv(4) * transpose(conjg(-dGR % matrix)) * nambuv(4)
  end function

  pure function propagator_keldysh(this) result(GK)
    !! Calculates the 4×4 Keldysh propagator G^K.
    class(propagator), intent(in) :: this     !! Propagator object
    type(nambu)                   :: GK       !! Propagator matrix
    type(nambu)                   :: GR, GA, H

    ! Calculate equilibrium propagators and the distribution
    GR = this % retarded()
    GA = this % advanced()
    H  = this % distribution()

    ! Use this to calculate the nonequilibrium propagator
    GK = GR*H - H*GA
  end function

  pure function propagator_keldysh_gradient(this, gauge) result(dGK)
    !! Calculates the 4×4 Keldysh propagator gradient ∇G^K. If an optional
    !! gauge field is specified, it returns the gauge-covariant gradient.
    class(propagator),     intent(in) :: this     !! Propagator object
    type(nambu), optional, intent(in) :: gauge    !! Optional gauge field
    type(nambu)                       :: dGK      !! Propagator gradient
    type(nambu)                       :: GR,  GA,  H
    type(nambu)                       :: dGR, dGA, dH

    ! Calculate equilibrium propagators and the distribution
    GR  = this % retarded()
    GA  = this % advanced()
    H   = this % distribution()

    ! Calculate the gradients of the matrix functions above
    dGR = this % retarded_gradient(gauge)
    dGA = this % advanced_gradient(gauge)
    dH  = this % distribution_gradient(gauge)

    ! Use this to calculate the nonequilibrium propagator gradient
    dGK = (dGR*H - H*dGA) + (GR*dH - dH*GA)
  end function

  pure function propagator_distribution(this) result(H)
    !! Calculates the 4×4 distribution function matrix H.
    class(propagator), intent(in) :: this   !! Propagator object
    type(nambu)                   :: H      !! Distribution matrix
    integer                       :: i

    ! Construct the distribution matrix from its Pauli-decomposition
    do i=0,7
      H = H + nambuv(i) * this % h(i)
    end do
  end function

  pure function propagator_distribution_gradient(this, gauge) result(dH)
    !! Calculates the 4×4 distribution function gradient ∇H. If an optional
    !! gauge field is specified, it returns the gauge-covariant gradient.
    class(propagator),     intent(in) :: this    !! Propagator object
    type(nambu), optional, intent(in) :: gauge   !! Optional gauge field
    type(nambu)                       :: dH      !! Distribution gradient
    integer                           :: i

    ! Construct the distribution matrix from its Pauli-decomposition
    do i=0,7
      dH = dH + nambuv(i) * this % dh(i)
    end do

    ! Construct the gauge-covariant terms
    if (present(gauge)) then
      associate (A  => gauge, H => this % distribution(), i  => (0.0_wp,1.0_wp))
        dH = dH - i*(A*H-H*A)
      end associate
    end if
  end function

  pure subroutine propagator_import_rvector(a, b)
    ! Defines assignment from a real vector to a propagator object.
    type(propagator), intent(inout) :: a
    real(wp),         intent(in)    :: b(32)

    ! Unpack the vector elements
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

  pure function propagator_supercurrent(this, gauge) result(J)
    !! Calculates the spectral supercurrents in the junction. The result is returned in the
    !! form of an 8-vector containing the charge, spin, heat, and spin-heat currents.
    class(propagator),     intent(in) :: this     !! Propagator object
    type(nambu), optional, intent(in) :: gauge    !! Optional gauge field
    real(wp), dimension(0:7)          :: J        !! Spectral supercurrent
    type(nambu)                       :: I        !! Matrix supercurrent
    type(nambu)                       :: H        !! Distribution function
    type(nambu)                       :: GR, dGR  !! Retarded propagator
    type(nambu)                       :: GA, dGA  !! Advanced propagator
    integer                           :: n

    ! Calculate the propagators
    H   = this % distribution()
    GR  = this % retarded()
    GA  = this % advanced()
    dGR = this % retarded_gradient(gauge)
    dGA = this % advanced_gradient(gauge)

    ! Calculate the matrix current
    I = (GR*dGR)*H - H*(GA*dGA)

    ! Pauli-decompose the current
    do n=0,7
      J(n) = re(trace((nambuv(4) * nambuv(n)) * I))/4
    end do
  end function

  pure function propagator_lossycurrent(this, gauge) result(J)
    !! Calculates the spectral dissipative currents in the junction. The result is returned in
    !! the form of an 8-vector containing the charge, spin, heat, and spin-heat currents.
    class(propagator),     intent(in) :: this     !! Propagator object
    type(nambu), optional, intent(in) :: gauge    !! Optional gauge field
    real(wp), dimension(0:7)          :: J        !! Spectral dissipative current
    type(nambu)                       :: I        !! Matrix dissipative current
    type(nambu)                       :: dH       !! Distribution function
    type(nambu)                       :: GR       !! Retarded propagator
    type(nambu)                       :: GA       !! Advanced propagator
    integer                           :: n

    ! Calculate the propagators
    dH = this % distribution_gradient(gauge)
    GR = this % retarded()
    GA = this % advanced()

    ! Calculate the matrix current
    I = dH - GR*dH*GA

    ! Pauli-decompose the current
    do n=0,7
      J(n) = re(trace((nambuv(4) * nambuv(n)) * I))/4
    end do
  end function

  pure function propagator_accumulation(this) result(Q)
    !! Calculates the spectral accumulations in the junction. The result is returned in the
    !! form of an 8-vector containing the charge, spin, heat, and spin-heat accumulations.
    class(propagator), intent(in) :: this     !! Propagator object
    real(wp), dimension(0:7)      :: Q        !! Spectral accumulation
    type(nambu)                   :: GK       !! Keldysh propagator
    integer                       :: n

    ! Calculate the propagator
    GK = this % keldysh()

    ! Pauli-decompose it
    do n=0,7
      Q(n) = -re(trace(nambuv(n) * GK))/4
    end do
  end function

  pure function propagator_correlation(this) result(r)
    !! Calculates the spectral pair-correlation function. This is useful e.g. for
    !! self-consistently calculating the superconducting gap in a superconductor.
    class(propagator), intent(in) :: this    !! Propagator object
    complex(wp)                   :: r       !! Spectral correlation
    type(nambu)                   :: GK      !! Keldysh propagator
    type(spin)                    :: f, ft   !! Anomalous propagators

    ! Calculate the propagator
    GK = this % keldysh()

    ! Extract the anomalous components
    f  = GK % matrix(1:2,3:4) 
    ft = GK % matrix(3:4,1:2)

    ! Trace out the singlet component
    r = trace((0.0_wp,-1.0_wp) * pauli2 * (f+conjg(ft)))/8
  end function

  pure function propagator_density(this) result(D)
    !! Calculates the spin-resolved local density of states.
    class(propagator), intent(in) :: this
    real(wp), dimension(0:7)      :: D
    type(nambu)                   :: GR
    type(spin)                    :: g, gt
 
    ! Extract the normal retarded propagator
    GR = this % retarded()
    g  = +GR % matrix(1:2,1:2) 
    gt = -GR % matrix(3:4,3:4) 

    ! Calculate the spin-resolved density of states,
    ! for both positive and negative energy values
    D(0:3) = re(trace(pauli * g))/2
    D(4:7) = re(trace(pauli * gt))/2
  end function

  pure subroutine propagator_decompose(this, f, ft, df, dft)
    !! Performs a singlet/triplet decomposition of the anomalous retarded propagators (f, f~)
    !! and their derivatives (∇f, ∇f~). Each of these are returned as a complex 4-vector,
    !! where f(0) corresponds to the singlet component and f(1:3) to the triplet component.
    !!
    !! @NOTE: This function can only be used safely in equilibrium systems.
    class(propagator),                     intent(in)  :: this  !! Propagator object
    complex(wp), dimension(0:3), optional, intent(out) :: f     !! Anomalous propagator (f )
    complex(wp), dimension(0:3), optional, intent(out) :: ft    !! Anomalous propagator (f~)
    complex(wp), dimension(0:3), optional, intent(out) :: df    !! Anomalous propagator gradient (∇f )
    complex(wp), dimension(0:3), optional, intent(out) :: dft   !! Anomalous propagator gradient (∇f~)
    type(spin),  dimension(0:3)                        :: P, Pt !  Pauli vectors used for the calculation
    type(nambu)                                        :: G     !  4×4 retarded propagator (G^R)

    ! Calculate the Pauli matrices necessary to invert f=[f·σ]iσ₂
    P  = ((0.0_wp,-0.5_wp) * pauli2) * pauli
    Pt = conjg(P)


    ! Calculate and decompose the propagator
    if (present(f) .or. present(ft)) then
      ! Extract the retarded propagator
      G  = this % retarded()

      ! Perform the singlet/triplet decomposition
      if (present(f )) then
        f   = trace(spin(+G % matrix(1:2, 3:4)) * P )
      end if
      if (present(ft)) then
        ft  = trace(spin(-G % matrix(3:4, 1:2)) * Pt)
      end if
    end if

    ! Calculate and decompose the propagator gradient
    if (present(df) .or. present(dft)) then
      ! Extract the retarded propagator gradient
      G  = this % retarded_gradient()

      ! Perform the singlet/triplet decomposition
      if (present(df )) then
        df   = trace(spin(+G % matrix(1:2, 3:4)) * P )
      end if
      if (present(dft)) then
        dft  = trace(spin(-G % matrix(3:4, 1:2)) * Pt)
      end if
    end if
  end subroutine
end module
