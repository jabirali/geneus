! This file defines a module containing 
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-10
! Updated: 2015-07-10

program mytest
  use module_precision
  use module_spin
  !use module_structures
  implicit none

  real(dp) :: a = 2
  type(spin) :: test

  print *,test
  test = spin([ 1., 2., 3., 4., 5., 6., 7., 8.])
  print *,test
  test = spin([ (2,1), (4,3), (6,5), (8,7) ])
  print *,test
  test = spin(reshape([ (1,2), (3,4), (5,6), (7,8) ], [2,2]))
  print *,test
  test = 2. * test * 2.
  print *,test
  test = (0,1) * test * (0,1)
  print *,test
  test = reshape([ (2,0), (0,0), (0,0), (2,0) ], [2,2]) * test * reshape([ (2,0), (0,0), (0,0), (2,0) ], [2,2])
  print *,test
end program
