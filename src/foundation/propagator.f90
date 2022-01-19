!> Author:   Jabir Ali Ouassou
!> Category: Foundation
!>
!> This module defines the type 'propagator', which represents a propagator
!> (i.e. Green's function) at a given position and energy. The equilibrium
!> parts (retarded and advanced) are represented via the Riccati parameters
!> γ, γ~ and their derivatives, while the nonequilibrium part (Keldysh) is
!> represented as traces of the distribution function and its derivatives.
!> These are together sufficient to reconstruct the full 8×8 propagator and
!> its derivatives, and can be used to calculate associated physical quantities
!> such as the density of states, charge currents, spin currents, heat currents,
!> spin-heat currents, and various nonequilibrium quasiparticle accumulations.

module propagator_m
    use :: math_m
    use :: spin_m
    use :: nambu_m
    private

    ! Public interface
    public propagator

    ! Type declaration
    type propagator
        ! Riccati parametrization of equilibrium propagators.
        type(spin) :: g    !! Riccati parameter γ
        type(spin) :: gt   !! Riccati parameter γ~
        type(spin) :: dg   !! Riccati parameter ∇γ
        type(spin) :: dgt  !! Riccati parameter ∇γ~
        type(spin) :: d2g  !! Riccati parameter ∇²γ
        type(spin) :: d2gt !! Riccati parameter ∇²γ~
        type(spin) :: N    !! Riccati normalization N
        type(spin) :: Nt   !! Riccati normalization N~

        ! Distribution-trace parametrization of nonequilibrium propagators.
        real(wp), dimension(0:7) :: h   = [1, 0, 0, 0, 0, 0, 0, 0] !! Distribution trace H
        real(wp), dimension(0:7) :: dh  = [0, 0, 0, 0, 0, 0, 0, 0] !! Distribution trace ∇H
        real(wp), dimension(0:7) :: d2h = [0, 0, 0, 0, 0, 0, 0, 0] !! Distribution trace ∇²H
    contains
        ! Accessors for the propagator matrices represented by this object
        procedure :: retarded             => propagator_retarded              !! Retarded propagator Gᴿ
        procedure :: retarded_gradient    => propagator_retarded_gradient     !! Retarded propagator ∇Gᴿ
        procedure :: retarded_laplacian   => propagator_retarded_laplacian    !! Retarded propagator ∇²Gᴿ

        procedure :: advanced             => propagator_advanced              !! Advanced propagator Gᴬ
        procedure :: advanced_gradient    => propagator_advanced_gradient     !! Advanced propagator ∇Gᴬ
        procedure :: advanced_laplacian   => propagator_advanced_laplacian    !! Advanced propagator ∇²Gᴬ

        procedure :: keldysh              => propagator_keldysh               !! Keldysh propagator Gᴷ
        procedure :: keldysh_gradient     => propagator_keldysh_gradient      !! Keldysh propagator ∇Gᴷ

        procedure :: distribution         => propagator_distribution          !! Distribution matrix H
        procedure :: distribution_gradien => propagator_distribution_gradient !! Distribution matrix ∇H

        ! Accessors for derived matrices used to solve the kinetic equations
        procedure :: dissipation          => propagator_dissipation           !! Dissipation matrix M
        procedure :: dissipation_gradient => propagator_dissipation_gradient  !! Dissipation matrix ∇M

        procedure :: condensate           => propagator_condensate            !! Condensate matrix Q
        procedure :: condensate_gradient  => propagator_condensate_gradient   !! Condensate matrix ∇Q

        procedure :: selfenergy1          => propagator_selfenergy1           !! Selfenergy matrix R₁
        procedure :: selfenergy2          => propagator_selfenergy2           !! Selfenergy matrix R₂

        ! Accessors for physical quantities that derive from the propagators
        procedure :: supercurrent         => propagator_supercurrent          !! Spectral supercurrents
        procedure :: lossycurrent         => propagator_lossycurrent          !! Spectral dissipative currents
        procedure :: accumulation         => propagator_accumulation          !! Spectral accumulations
        procedure :: correlation          => propagator_correlation           !! Spectral correlations
        procedure :: density              => propagator_density               !! Spin-resolved density of states

        ! Miscellaneous utiliy functions for working with propagator objects
        procedure :: save                 => propagator_save                  !! Export Riccati parameters
        procedure :: load                 => propagator_load                  !! Import Riccati parameters
    end type

    ! Type constructor
    interface propagator
        module procedure propagator_construct_vacuum, propagator_construct_riccati, propagator_construct_bcs
    end interface
contains
    pure function propagator_construct_vacuum() result(this)
    !!  Construct a vacuum propagator, i.e. a propagator which satisfies G=0.
        type(propagator) :: this !! Constructed object

        continue
    end function

    pure function propagator_construct_riccati(g, gt, dg, dgt) result(this)
    !!  Construct an arbitrary state by explicitly providing Riccati parameters.
    !!  Unspecified Riccati parameters default to zero due to spin constructors.
    !!  The distribution function defaults to equilibrium at zero temperature.
        type(propagator)                 :: this !! Constructed object
        type(spin),           intent(in) :: g    !! Riccati parameter γ
        type(spin),           intent(in) :: gt   !! Riccati parameter γ~
        type(spin), optional, intent(in) :: dg   !! Riccati parameter ∇γ
        type(spin), optional, intent(in) :: dgt  !! Riccati parameter ∇γ~

        ! Copy Riccati parameters into the new object
        this%g  = g
        this%gt = gt

        ! Copy Riccati derivatives into the new object
        if (present(dg))  this%dg  = dg
        if (present(dgt)) this%dgt = dgt

        ! Update the normalization matrices
        associate (g => this%g, gt => this%gt, N => this%N, Nt => this%Nt)
            N  = inverse(pauli0 - g*gt)
            Nt = inverse(pauli0 - gt*g)
        end associate
    end function

    pure function propagator_construct_bcs(energy, gap) result(this)
    !!  Constructs the state of a a BCS superconductor at a given energy, which
    !!  may have an imaginary term representing inelastic scattering. The second
    !!  argument 'gap' is used to provide the superconducting order parameter Δ.
    !!  The distribution function defaults to equilibrium at zero temperature.
        type(propagator)        :: this   !! Constructed object
        complex(wp), intent(in) :: energy !! Quasiparticle energy
        complex(wp), intent(in) :: gap    !! Order parameter

        real(wp)    :: p
        complex(wp) :: t, u
        complex(wp) :: a, b

        ! Calculate the superconducting gap and phase
        u = abs(gap)/energy
        p = arg(gap)

        ! Calculate the θ-parameter
        t = (log(1 + u) - log(1 - u))/2

        ! Calculate the scalar Riccati parameters a and b
        a =  (exp(+t) - exp(-t))/(2 + exp(+t) + exp(-t))*exp((0, +1)*p)
        b = -(exp(+t) - exp(-t))/(2 + exp(+t) + exp(-t))*exp((0, -1)*p)

        ! Calculate the matrix Riccati parameters γ and γ~
        this%g  = a*((0.0_wp, 1.0_wp)*pauli2)
        this%gt = b*((0.0_wp, 1.0_wp)*pauli2)

        ! Update the normalization matrices
        this%N  = inverse(pauli0 - this%g*this%gt)
        this%Nt = inverse(pauli0 - this%gt*this%g)
    end function

    pure function propagator_retarded(this) result(GR)
    !!  Calculates the 4×4 retarded propagator Gᴿ.
        class(propagator), intent(in) :: this !! Propagator object
        type(nambu)                   :: GR   !! Retarded propagator

        ! Construct the propagator from the Riccati parameters
        associate (g => this%g, gt => this%gt, &
                   N => this%N, Nt => this%Nt, &
                   I => pauli0, M  => GR%matrix)
            M(1:2, 1:2) = (+1.0_wp)*N*(I + g*gt)
            M(1:2, 3:4) = (+2.0_wp)*N*g
            M(3:4, 1:2) = (-2.0_wp)*Nt*gt
            M(3:4, 3:4) = (-1.0_wp)*Nt*(I + gt*g)
        end associate
    end function

    pure function propagator_retarded_gradient(this, gauge) result(dGR)
    !!  Calculates the 4×4 retarded propagator gradient ∇Gᴿ. If an optional
    !!  gauge field is specified, it returns the gauge-covariant gradient.
        class(propagator),     intent(in) :: this  !! Propagator object
        type(nambu), optional, intent(in) :: gauge !! Optional gauge field
        type(nambu)                       :: dGR   !! Retarded propagator gradient

        ! Construct the propagator from the Riccati parameters
        associate (g  => this%g,  gt  => this%gt,  &
                   dg => this%dg, dgt => this%dgt, &
                   N  => this%N,  Nt  => this%Nt,  &
                   I  => pauli0,  M   => dGR%matrix)
            M(1:2, 1:2) = (+2.0_wp)*N*(dg*gt + g*dgt)*N
            M(1:2, 3:4) = (+2.0_wp)*N*(dg + g*dgt*g)*Nt
            M(3:4, 1:2) = (-2.0_wp)*Nt*(dgt + gt*dg*gt)*N
            M(3:4, 3:4) = (-2.0_wp)*Nt*(dgt*g + gt*dg)*Nt
        end associate

        ! Construct the gauge-covariant terms
        if (present(gauge)) then
            associate (A => gauge, GR => this%retarded(), i => (0.0_wp, 1.0_wp))
                dGR = dGR - i*(A*GR - GR*A)
            end associate
        end if
    end function

    pure function propagator_retarded_laplacian(this) result(d2GR)
    !!  Calculates the 4×4 retarded propagator gradient ∇²Gᴿ.
    !!
    !!  @TODO:
    !!    Implement support for gauge-covariant laplacians.
        class(propagator), intent(in) :: this !! Propagator object
        type(nambu)                   :: d2GR !! Retarded propagator laplacian

        type(spin) :: D, Dt, dD, dDt, F, Ft, dF, dFt

        ! Construct the propagator from the Riccati parameters
        associate (g   => this%g,   gt   => this%gt,   &
                   dg  => this%dg,  dgt  => this%dgt,  &
                   d2g => this%d2g, d2gt => this%d2gt, &
                   N   => this%N,   Nt   => this%Nt,   &
                   I   => pauli0,   M    => d2GR%matrix)

            ! Calculate 1st-derivative auxiliary matrices
            D  = dg*gt + g*dgt
            Dt = dgt*g + gt*dg
            F  = dg  + g*dgt*g
            Ft = dgt + gt*dg*gt

            ! Calculate 2nd-derivative auxiliary matrices
            dD  = d2g*gt + g*d2gt + 2.0_wp*dg*dgt
            dDt = d2gt*g + gt*d2g + 2.0_wp*dgt*dg
            dF  = d2g  + g*d2gt*g  + dg*dgt*g  + g*dgt*dg
            dFt = d2gt + gt*d2g*gt + dgt*dg*gt + gt*dg*dgt

            ! Calculate the propagator matrix
            M(1:2, 1:2) = (+2.0_wp)*N*(dD + 2.0_wp*D*N*D)*N
            M(1:2, 3:4) = (+2.0_wp)*N*(dF + D*N*F + F*Nt*Dt)*Nt
            M(3:4, 1:2) = (-2.0_wp)*Nt*(dFt + Dt*Nt*Ft + Ft*N*D)*N
            M(3:4, 3:4) = (-2.0_wp)*Nt*(dDt + 2.0_wp*Dt*Nt*Dt)*Nt
        end associate
    end function

    pure function propagator_advanced(this) result(GA)
    !!  Calculates the 4×4 advanced propagator Gᴬ.
        class(propagator), intent(in) :: this !! Propagator object
        type(nambu)                   :: GA   !! Advanced propagator

        type(nambu) :: GR

        ! Calculate the retarded propagator
        GR = this%retarded()

        ! Use the identity GA = -τ₃GR†τ₃
        GA = nambuv(4)*transpose(conjg(-GR%matrix))*nambuv(4)
    end function

    pure function propagator_advanced_gradient(this, gauge) result(dGA)
    !!  Calculates the 4×4 advanced propagator gradient ∇Gᴬ. If an optional
    !!  gauge field is specified, it returns the gauge-covariant gradient.
        class(propagator),     intent(in) :: this  !! Propagator object
        type(nambu), optional, intent(in) :: gauge !! Optional gauge field
        type(nambu)                       :: dGA   !! Advanced propagator gradient

        type(nambu) :: dGR

        ! Calculate the retarded propagator gradient
        dGR = this%retarded_gradient(gauge)

        ! Use the identity GA = -τ₃GR†τ₃
        dGA = nambuv(4)*transpose(conjg(-dGR%matrix))*nambuv(4)
    end function

    pure function propagator_advanced_laplacian(this) result(d2GA)
    !!  Calculates the 4×4 retarded propagator gradient ∇²Gᴬ.
    !!
    !!  @TODO:
    !!    Implement support for gauge-covariant laplacians.
        class(propagator), intent(in) :: this !! Propagator object
        type(nambu)                   :: d2GA !! Advanced propagator laplacian

        type(nambu) :: d2GR

        ! Calculate the retarded propagator laplacian
        d2GR = this%retarded_laplacian()

        ! Use the identity GA = -τ₃GR†τ₃
        d2GA = nambuv(4)*transpose(conjg(-d2GR%matrix))*nambuv(4)
    end function

    pure function propagator_keldysh(this) result(GK)
    !!  Calculates the 4×4 Keldysh propagator Gᴷ.
        class(propagator), intent(in) :: this !! Propagator object
        type(nambu)                   :: GK   !! Propagator matrix

        type(nambu) :: GR, GA, H

        ! Calculate equilibrium propagators and the distribution
        GR = this%retarded()
        GA = this%advanced()
        H  = this%distribution()

        ! Use this to calculate the nonequilibrium propagator
        GK = GR*H - H*GA
    end function

    pure function propagator_keldysh_gradient(this, gauge) result(dGK)
    !!  Calculates the 4×4 Keldysh propagator gradient ∇Gᴷ. If an optional
    !!  gauge field is specified, it returns the gauge-covariant gradient.
        class(propagator),     intent(in) :: this  !! Propagator object
        type(nambu), optional, intent(in) :: gauge !! Optional gauge field
        type(nambu)                       :: dGK   !! Propagator gradient

        type(nambu) :: GR, GA, H
        type(nambu) :: dGR, dGA, dH

        ! Calculate equilibrium propagators and the distribution
        GR = this%retarded()
        GA = this%advanced()
        H  = this%distribution()

        ! Calculate the gradients of the matrix functions above
        dGR = this%retarded_gradient(gauge)
        dGA = this%advanced_gradient(gauge)
        dH  = this%distribution_gradient(gauge)

        ! Use this to calculate the nonequilibrium propagator gradient
        dGK = (dGR*H - H*dGA) + (GR*dH - dH*GA)
    end function

    pure function propagator_distribution(this) result(H)
    !!  Calculates the 4×4 distribution function matrix H.
        class(propagator), intent(in) :: this !! Propagator object
        type(nambu)                   :: H    !! Distribution matrix

        integer :: i

        ! Construct the distribution matrix from its Pauli-decomposition
        do i = 0, 7
            H = H + nambuv(i)*this%h(i)
        end do
    end function

    pure function propagator_distribution_gradient(this, gauge) result(dH)
    !!  Calculates the 4×4 distribution function gradient ∇H. If an optional
    !!  gauge field is specified, it returns the gauge-covariant gradient.
        class(propagator),     intent(in) :: this  !! Propagator object
        type(nambu), optional, intent(in) :: gauge !! Optional gauge field
        type(nambu)                       :: dH    !! Distribution gradient

        integer :: i

        ! Construct the distribution matrix from its Pauli-decomposition
        do i = 0, 7
            dH = dH + nambuv(i)*this%dh(i)
        end do

        ! Construct the gauge-covariant terms
        if (present(gauge)) then
            associate (A => gauge, H => this%distribution(), i => (0.0_wp, 1.0_wp))
                dH = dH - i*(A*H - H*A)
            end associate
        end if
    end function

    pure function propagator_supercurrent(this, gauge) result(J)
    !!  Calculates the spectral supercurrents in the junction. The result is an
    !!  8-vector encoding respectively charge, spin, heat, spin-heat currents.
        class(propagator),     intent(in) :: this  !! Propagator object
        type(nambu), optional, intent(in) :: gauge !! Optional gauge field
        real(wp), dimension(0:7)          :: J     !! Spectral supercurrent

        type(nambu) :: I, H, GR, dGR, GA, dGA
        integer     :: n

        ! Calculate the propagators
        H = this%distribution()
        GR = this%retarded()
        GA = this%advanced()
        dGR = this%retarded_gradient(gauge)
        dGA = this%advanced_gradient(gauge)

        ! Calculate the matrix current
        I = (GR*dGR)*H - H*(GA*dGA)

        ! Pauli-decompose the current
        do n = 0, 7
            J(n) = re(trace((nambuv(4)*nambuv(n))*I))/8
        end do
    end function

    pure function propagator_lossycurrent(this, gauge) result(J)
    !!  Calculates the spectral dissipative currents in the junction. The result
    !!  is an 8-vector containing the charge, spin, heat, and spin-heat currents.
        class(propagator),     intent(in) :: this  !! Propagator object
        type(nambu), optional, intent(in) :: gauge !! Optional gauge field
        real(wp), dimension(0:7)          :: J     !! Spectral dissipative current

        type(nambu) :: I, dH, GR, GA
        integer     :: n

        ! Calculate the propagators
        dH = this%distribution_gradient(gauge)
        GR = this%retarded()
        GA = this%advanced()

        ! Calculate the matrix current
        I = dH - GR*dH*GA

        ! Pauli-decompose the current
        do n = 0, 7
            J(n) = re(trace((nambuv(4)*nambuv(n))*I))/8
        end do
    end function

    pure function propagator_accumulation(this) result(Q)
    !!  Calculates the spectral accumulations. The result is an 8-vector
    !!  containing the charge, spin, heat, and spin-heat accumulations.
        class(propagator), intent(in) :: this !! Propagator object
        real(wp), dimension(0:7)      :: Q    !! Spectral accumulation

        type(nambu) :: GK
        integer     :: n

        ! Calculate the propagator
        GK = this%keldysh()

        ! Pauli-decompose it
        do n = 0, 7
            Q(n) = -re(trace(nambuv(n)*GK))/8
        end do
    end function

    pure function propagator_correlation(this) result(r)
    !!  Calculates the spectral pair-correlation function. This is useful when
    !!  self-consistently calculating the superconducting gap in a superconductor.
        class(propagator), intent(in) :: this !! Propagator object
        complex(wp)                   :: r    !! Spectral correlation

        type(nambu) :: GK
        type(spin)  :: f, ft

        ! Calculate the propagator
        GK = this%keldysh()

        ! Extract the anomalous components
        f  = GK%matrix(1:2, 3:4)
        ft = GK%matrix(3:4, 1:2)

        ! Trace out the singlet component
        r = trace((0.0_wp, -1.0_wp)*pauli2*(f + conjg(ft)))/8
    end function

    pure function propagator_density(this) result(D)
    !!  Calculates the spin-resolved local density of states.
        class(propagator), intent(in) :: this
        real(wp), dimension(0:7)      :: D

        type(nambu) :: GR
        type(spin)  :: g, gt

        ! Extract the normal retarded propagator
        GR = this%retarded()
        g  = +GR%matrix(1:2, 1:2)
        gt = -GR%matrix(3:4, 3:4)

        ! Calculate the spin-resolved density of states,
        ! for both positive and negative energy values
        D(0:3) = re(trace(pauli*g))/2
        D(4:7) = re(trace(pauli*gt))/2
    end function

    pure elemental subroutine propagator_save(this, other)
    !!  Defines a function for exporting Riccati parameters.
        class(propagator), intent(inout) :: this
        class(propagator), intent(inout) :: other

        ! Copy all the Riccati parameters
        other%g    = this%g
        other%gt   = this%gt
        other%dg   = this%dg
        other%dgt  = this%dgt
        other%d2g  = this%d2g
        other%d2gt = this%d2gt
        other%N    = this%N
        other%Nt   = this%Nt
    end subroutine

    pure elemental subroutine propagator_load(this, other)
    !!  Defines a function for importing Riccati parameters.
        class(propagator), intent(inout) :: this
        class(propagator), intent(inout) :: other

        ! Copy all the Riccati parameters
        this%g    = other%g
        this%gt   = other%gt
        this%dg   = other%dg
        this%dgt  = other%dgt
        this%d2g  = other%d2g
        this%d2gt = other%d2gt
        this%N    = other%N
        this%Nt   = other%Nt
    end subroutine

    pure function propagator_dissipation(this) result(M)
    !!  Calculates the dissipation matrix M = ∂J/∂H', where J is the
    !!  current and H' is the gradient of the distribution function.
        class(propagator), intent(in)    :: this
        complex(wp), dimension(0:7, 0:7) :: M

        type(nambu), dimension(0:7) :: N
        type(nambu)                 :: GR, GA
        integer                     :: i, j

        ! Memoize the basis matrices
        do i = 0, 7
            N(i) = nambuv(i)
        end do

        ! Construct the propagator matrices
        GR = this%retarded()
        GA = this%advanced()

        ! Construct the dissipation matrix
        do i = 0, 7
            do j = 0, 7
                M(i, j) = trace(N(i)*N(j) - N(i)*GR*N(j)*GA)/4
            end do
        end do
    end function

    pure function propagator_dissipation_gradient(this) result(dM)
    !!  Calculates the gradient of the dissipation matrix M'.
        class(propagator), intent(in)    :: this
        complex(wp), dimension(0:7, 0:7) :: dM

        type(nambu), dimension(0:7) :: N
        type(nambu)                 :: GR, GA, dGR, dGA
        integer                     :: i, j

        ! Memoize the basis matrices
        do i = 0, 7
            N(i) = nambuv(i)
        end do

        ! Construct the propagator matrices
        GR = this%retarded()
        GA = this%advanced()

        ! Construct the propagator gradients
        dGR = this%retarded_gradient()
        dGA = this%advanced_gradient()

        ! Construct the dissipation matrix
        do j = 0, 7
            do i = 0, 7
                dM(i, j) = -trace(N(i)*dGR*N(j)*GA + N(i)*GR*N(j)*dGA)/4
            end do
        end do
    end function

    pure function propagator_condensate(this) result(Q)
    !!  Calculates the condensate matrix Q = ∂J/∂H, where J is the
    !!  current and H is the nonequilibrium distribution function.
        class(propagator), intent(in)    :: this
        complex(wp), dimension(0:7, 0:7) :: Q

        type(nambu), dimension(0:7) :: N
        type(nambu)                 :: GR, GA, dGR, dGA
        integer                     :: i, j

        ! Memoize the basis matrices
        do i = 0, 7
            N(i) = nambuv(i)
        end do

        ! Construct the propagator matrices
        GR = this%retarded()
        GA = this%advanced()

        ! Construct the propagator gradients
        dGR = this%retarded_gradient()
        dGA = this%advanced_gradient()

        ! Construct the condensate matrix
        do j = 0, 7
            do i = 0, 7
                Q(i, j) = trace(N(j)*N(i)*GR*dGR - N(i)*N(j)*GA*dGA)/4
            end do
        end do
    end function

    pure function propagator_condensate_gradient(this) result(dQ)
    !!  Calculates the gradient of the condensate matrix Q'.
        class(propagator), intent(in)    :: this
        complex(wp), dimension(0:7, 0:7) :: dQ

        type(nambu), dimension(0:7) :: N
        type(nambu)                 :: GR, GA, dGR, dGA, d2GR, d2GA
        integer                     :: i, j

        ! Memoize the basis matrices
        do i = 0, 7
            N(i) = nambuv(i)
        end do

        ! Construct the propagator matrices
        GR = this%retarded()
        GA = this%advanced()

        ! Construct the propagator gradients
        dGR = this%retarded_gradient()
        dGA = this%advanced_gradient()

        ! Construct the propagator laplacians
        d2GR = this%retarded_laplacian()
        d2GA = this%advanced_laplacian()

        ! Construct the condensate matrix
        do j = 0, 7
            do i = 0, 7
                dQ(i, j) = trace(N(j)*N(i)*(dGR*dGR + GR*d2GR) - N(i)*N(j)*(dGA*dGA + GA*d2GA))/4
            end do
        end do
    end function

    pure function propagator_selfenergy1(this, S) result(R)
    !!  Calculates the 1st-order self-energy contribution to the kinetic equations.
        class(propagator), intent(in)    :: this
        type(nambu),       intent(in)    :: S
        complex(wp), dimension(0:7, 0:7) :: R

        type(nambu), dimension(0:7) :: N
        type(nambu)                 :: GR, GA
        integer                     :: i, j

        ! Memoize the basis matrices
        do i = 0, 7
            N(i) = nambuv(i)
        end do

        ! Construct the propagator matrices
        GR = this%retarded()
        GA = this%advanced()

        ! Construct the self-energy matrix
        do j = 0, 7
            do i = 0, 7
                R(i, j) = (0.00, 0.25)*trace((N(i)*S - S*N(i))*(GR*N(j) - N(j)*GA))
            end do
        end do
    end function

    pure function propagator_selfenergy2(this, S) result(R)
    !!  Calculates the 2nd-order self-energy contribution to the kinetic equations.
        class(propagator), intent(in)    :: this
        type(nambu),       intent(in)    :: S
        complex(wp), dimension(0:7, 0:7) :: R

        type(nambu), dimension(0:7) :: N
        type(nambu)                 :: GR, GA
        integer                     :: i, j

        ! Memoize the basis matrices
        do i = 0, 7
            N(i) = nambuv(i)
        end do

        ! Construct the propagator matrices
        GR = this%retarded()
        GA = this%advanced()

        ! Construct the self-energy matrix
        do j = 0, 7
            do i = 0, 7
                R(i, j) = (0.00, 0.25)*trace((N(i)*S - S*N(i))*(GR*S*GR*N(j) - N(j)*GA*S*GA + GR*(N(j)*S - S*N(j))*GA))
            end do
        end do
    end function
end module
