!> Author:   Jabir Ali Ouassou
!> Category: Materials
!>
!> This submodule is included by conductor.f, and contains the equations which model spin-active interfaces.
!>
!> @TODO: Rewrite using the new nambu.f library, replacing e.g. diag(m·σ,m·σ*) with m*nambuv(1:3).
!>        Also, we may then for brevity replace matmul(G,matmul(M,G)) with G*M*G, and so on.

module spinactive_m
  use :: propagator_m
  use :: material_m
  use :: math_m
  use :: spin_m
  use :: nambu_m

  ! Type declarations
  type, public :: spinactive
    real(wp)                 :: conductance   = 0.0      !! Interfacial conductance
    real(wp)                 :: polarization  = 0.0      !! Interfacial spin-polarization
    real(wp)                 :: spinmixing    = 0.0      !! Interfacial 1st-order spin-mixing
    real(wp)                 :: secondorder   = 0.0      !! Interfacial 2nd-order spin-mixing
    real(wp), dimension(1:3) :: magnetization = [0,0,1]  !! Interfacial magnetization direction
    real(wp), dimension(1:3) :: misalignment  = [0,0,0]  !! Interfacial magnetization misalignment
  
    type(nambu), private     :: M                        !! Magnetization matrix (transmission)
    type(nambu), private     :: M0                       !! Magnetization matrix (reflection, this  side)
    type(nambu), private     :: M1                       !! Magnetization matrix (reflection, other side)
  contains
    procedure :: diffusion_current    => spinactive_diffusion_current
    procedure :: update_prehook       => spinactive_update_prehook
  end type
contains
  pure subroutine spinactive_update_prehook(this)
    !! Updates the internal variables associated with spin-active interfaces.
    class(spinactive), intent(inout) :: this 
  
    ! Process transmission properties
    this % M = nambuv(this % magnetization)
  
    ! Default reflection properties match transmission properties
    this % M0 = this % M
    this % M1 = this % M
  
    ! Process reflection properties (this side)
    if (nonzero(this % misalignment)) then
      this % M0 = nambuv(this % misalignment)
    end if
  
    ! ! Process reflection properties (other side)
    ! if (associated(this % material_a)) then
    !   select type (other => this % material_a)
    !     class is (conductor)
    !       call update_magnetization(this % spinactive_a % M1, other % spinactive_b % misalignment)
    !   end select
    ! end if
  
    ! if (associated(this % material_b)) then
    !   select type (other => this % material_b)
    !     class is (conductor)
    !       call update_magnetization(this % spinactive_b % M1, other % spinactive_a % misalignment)
    !   end select
    ! end if
  end subroutine
  
  pure function spinactive_diffusion_current(this, G0, G1) result(I)
    !! Calculate the matrix current at an interface with spin-active properties. The equations
    !! implemented here should be valid for an arbitrary interface polarization, and up to 2nd
    !! order in the transmission probabilities and spin-mixing angles of the interface. 
    !! @TODO: Shortcut-evaluation for nonmagnetic interfaces.
    class(spinactive), intent(in) :: this  !! Bouundary object
    type(nambu),       intent(in) :: G0, G1      !! Propagator matrices
    type(nambu)                   :: S0, S1      !! Matrix expressions
    type(nambu)                   :: I           !! Matrix current
  
    ! Evaluate the 1st-order matrix functions
    S0 = spinactive_current1_transmission(G1)
    S1 = spinactive_current1_reflection()
  
    ! Evaluate the 1st-order matrix current
    associate(S => S0 + S1)
      I  = (G0*S - S*G0)
    end associate
  
    ! Calculate the 2nd-order contributions to the matrix current. Note that we make a
    ! number of simplifications in this implementation. In particular, we assume that 
    ! all interface parameters except the magnetization directions are equal on both
    ! sides of the interface. We also assume that the spin-mixing angles and tunneling
    ! probabilities of different channels have standard deviations that are much smaller
    ! than their mean values, which reduces the number of new fitting parameters to one.
  
    if (abs(this % secondorder) > 0) then
      ! Evaluate the 1st-order matrix functions
      associate(M1 => this % M1)
        S1 = spinactive_current1_transmission(G1*M1*G1 - M1) 
      end associate
  
      ! Evaluate the 2nd-order matrix current
      I = I                                     &
        + spinactive_current2_transmission()    &
        + spinactive_current2_crossterms()      &
        + spinactive_current2_reflection()
    end if

    ! Scale the answers based on conductance
    I = (this % conductance/2) * I
  contains
    pure function spinactive_current1_transmission(G) result(F)
      !! Calculate the 1st-order transmission terms in the matrix current commutator.
      type(nambu), intent(in) :: G
      type(nambu)             :: F
      real(wp) :: Pr, Pp, Pm
  
      associate(P => this % polarization, M => this % M)
        Pr = sqrt(1 - P**2)
        Pp = 1 + Pr
        Pm = 1 - Pr
  
        F = G + (P/Pp)*(M*G+G*M) + (Pm/Pp)*(M*G*M)
      end associate
    end function
  
    pure function spinactive_current1_reflection() result(F)
      !! Calculate the 1st-order spin-mixing terms in the matrix current commutator.
      type(nambu) :: F
  
      associate(Q => this % spinmixing, M0 => this % M0)
        F = ((0,-1)*Q) * M0
      end associate
    end function
  
    pure function spinactive_current2_transmission() result(I)
      !! Calculate the 2nd-order transmission terms in the matrix current.
      type(nambu) :: I
  
      associate(R => this % secondorder, Q => this % spinmixing)
        I = (0.50*R/Q) * (S0*G0*S0)
      end associate
    end function
  
    pure function spinactive_current2_reflection() result(I)
      !! Calculate the 2nd-order spin-mixing terms in the matrix current.
      type(nambu) :: I
  
      associate(R => this % secondorder, Q => this % spinmixing, M0 => this % M0)
        associate(U => M0*G0*M0)
          I = (0.25*R*Q) * (G0*U - U*G0)
        end associate
      end associate
    end function
  
    pure function spinactive_current2_crossterms() result(I)
      !! Calculate the 2nd-order cross-terms in the matrix current.
      type(nambu) :: I
  
      associate(R => this % secondorder, M0 => this % M0)
        associate(U => S0*G0*M0 + M0*G0*S0 + S1)
          I = ((0.00,0.25)*R) * (G0*U - U*G0)
        end associate
      end associate
    end function
  end function
end module
