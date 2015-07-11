! This module defines the data type 'state', which can be used to represent the physical state of a
! conductor at a given position and energy. This is done by storing the Riccati parameters γ and γ~,
! and their first derivatives dγ/dx and dγ~/dx, which is equivalent to knowledge about the retarded 
! Green's function at the given location and energy. To make it easier to interact with ODE solvers, 
! the assignment operator is overloaded in such a way that 'state' becomes isomorphic to real(32).
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-11
! Updated: 2015-07-11

module module_state
  use module_precision
  use module_spin
  implicit none

  ! Class declaration
  type state
    type(spin) :: g      ! Riccati parameter γ  
    type(spin) :: gt     ! Riccati parameter γ~
    type(spin) :: dg     ! Derivative dγ /dx
    type(spin) :: dgt    ! Derivative dγ~/dx
  end type

  ! Class constructor
  interface state
    module procedure state_construct_bcs
  end interface

  ! Assignment operator
  interface assignment(=)
    module procedure state_import_rvector, state_export_rvector
  end interface
contains
  pure function state_construct_bcs(energy, gap) result(this)
    ! Constructs a state corresponding to a BCS superconductor at some given energy, which may have an imaginary
    ! term representing inelastic scattering. The second argument 'gap' is the superconducting order parameter Δ.
    type(state)             :: this      ! State object that will be constructed
    complex(dp), intent(in) :: energy    ! Quasiparticle energy (including inelastic scattering contribution)
    complex(dp), intent(in) :: gap       ! Superconducting order parameter (including superconducting phase)

    ! Calculate the θ-parameter t, and then the scalar Riccati parameters a and b
    complex(dp)             :: t, a, b 
    t = atanh(gap/energy)
    a = sinh(t)/(1+cosh(t))
    b = conjg(a)

    ! Fill out the matrix Riccati parameters g and gt
    this%g  = [(0.0_dp,0.0_dp),  a, -a, (0.0_dp,0.0_dp)]
    this%gt = [(0.0_dp,0.0_dp), -b,  b, (0.0_dp,0.0_dp)]
  end function

  pure subroutine state_export_rvector(a, b)
    ! Defines assignment from a state object to a real vector
    real(dp),    intent(out) :: a(32)
    type(state), intent(in)  :: b

    a( 1: 8) = b%g
    a( 9:16) = b%gt
    a(17:24) = b%dg
    a(25:32) = b%dgt
  end subroutine

  pure subroutine state_import_rvector(a, b)
    ! Defines assignment from a real vector to a state object
    type(state), intent(out) :: a
    real(dp),    intent(in)  :: b(32)

    a%g   = b( 1: 8) 
    a%gt  = b( 9:16) 
    a%dg  = b(17:24) 
    a%dgt = b(25:32) 
  end subroutine
end module

program test
  use module_precision
  use module_spin
  use module_state
  type(state) :: mystate, newstate
  real(dp)    :: vec(32)

  print *,mystate

  mystate = state(spin(1.0_dp),spin(2.0_dp),spin(3.0_dp),spin(4.0_dp))

  print *,mystate%g
  print *,mystate%gt
  print *,mystate%dg
  print *,mystate%dgt
  
  print *,'Result:'
  vec = mystate
  print *,vec
  newstate = mystate
  print *,newstate
  
  print *,'New:'
  mystate = state( (1.2_dp,0.001_dp), (1.0_dp,0.0_dp) )
  print *,mystate
end program
