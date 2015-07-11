module module_state
  use module_precision
  use module_spin
  implicit none

  ! Class declaration
  type state
    type(spin) :: g   = state(0_dp) ! Riccati parameter γ  
    type(spin) :: gt  = state(0_dp) ! Riccati parameter γ~
    type(spin) :: dg  = state(0_dp) ! Derivative dγ /dx
    type(spin) :: dgt = state(0_dp) ! Derivative dγ~/dx
  end type

  ! Class constructor
  !interface state
  !  module procedure state_construct_superconductor, 
  !end interface

  ! writew

  ! Assignment operator
  interface assignment(=)
    module procedure state_import_rvector, state_export_rvector
  end interface
contains
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
end program
