! This program is used to test the data structures and procedures declared in 'modules/*.f90'.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-13
! Updated: 2015-07-13

program test_core
  use module_assert
  use module_spin

  real(dp)    :: rs, rv(8)
  complex(dp) :: cs, cv(4), cm(2,2)
  type(spin)  :: p, q, r


  ! Calibrate the test procedures in 'module_assert'
  call section('Calibration of assert subroutines')
  call subsection('Subroutines operating on a single value:')
  call assert(.true.,                'Logical true                                     SUCCESS == ')
  call assert(.false.,               'Logical false                                    FAILURE == ')
  call assert(0,                     'Integer 0                                        SUCCESS == ')
  call assert(1,                     'Integer 1                                        FAILURE == ')
  call assert(0.0_dp,                'Real    0                                        SUCCESS == ')
  call assert(1.0_dp,                'Real    1                                        FAILURE == ')
  call assert((0.0_dp,0.0_dp),       'Complex 0                                        SUCCESS == ')
  call assert((0.0_dp,1.0_dp),       'Complex i                                        FAILURE == ')
  call assert(spin(0.0_dp),          'Spin matrix 0                                    SUCCESS == ')
  call assert(spin(1.0_dp),          'Spin matrix 1                                    FAILURE == ')
  call assert(spin((0.0_dp,1.0_dp)), 'Spin matrix i                                    FAILURE == ')

  call subsection('Subroutines operating on arrays:')
  call assert([.true., .true.],                  'Logical [true, true ]                            SUCCESS == ')
  call assert([.true., .false.],                 'Logical [true, false]                            FAILURE == ')
  call assert([.false.,.false.],                 'Logical [false,false]                            FAILURE == ')
  call assert([0, 0],                            'Integer [0,0]                                    SUCCESS == ')
  call assert([1, 0],                            'Integer [1,0]                                    FAILURE == ')
  call assert([1, 1],                            'Integer [1,1]                                    FAILURE == ')
  call assert([0.0_dp,0.0_dp],                   'Real    [0,0]                                    SUCCESS == ')
  call assert([1.0_dp,0.0_dp],                   'Real    [1,0]                                    FAILURE == ')
  call assert([1.0_dp,1.0_dp],                   'Real    [1,1]                                    FAILURE == ')
  call assert([(0.0_dp,0.0_dp),(0.0_dp,0.0_dp)], 'Complex [0,0]                                    SUCCESS == ')
  call assert([(0.0_dp,1.0_dp),(0.0_dp,0.0_dp)], 'Complex [i,0]                                    FAILURE == ')
  call assert([(0.0_dp,1.0_dp),(0.0_dp,1.0_dp)], 'Complex [i,i]                                    FAILURE == ')



  ! Test the properties of the Pauli matrices in 'module_spin'
  call section('Pauli matrices')
  call subsection('Testing output:')
  call pauli0%print(' pauli0')
  call pauli1%print(' pauli1')
  call pauli2%print(' pauli2')
  call pauli3%print(' pauli3')

  call subsection('Testing squares:')
  call assert(pauli0*pauli0 - pauli0, 'pauli0^2 = pauli0')
  call assert(pauli1*pauli1 - pauli0, 'pauli1^2 = pauli0')
  call assert(pauli2*pauli2 - pauli0, 'pauli2^2 = pauli0')
  call assert(pauli3*pauli3 - pauli0, 'pauli3^2 = pauli0')

  call subsection('Testing inverses:')
  call assert(pauli0%inv() - pauli0, 'inv(pauli0) = pauli0')
  call assert(pauli1%inv() - pauli1, 'inv(pauli1) = pauli1')
  call assert(pauli2%inv() - pauli2, 'inv(pauli2) = pauli2')
  call assert(pauli3%inv() - pauli3, 'inv(pauli3) = pauli3')

  call subsection('Testing traces:')
  call assert(pauli0%trace() - (2.0_dp,0.0_dp),  'trace(pauli0) = 2')
  call assert(pauli1%trace() - (0.0_dp,0.0_dp),  'trace(pauli1) = 0')
  call assert(pauli2%trace() - (0.0_dp,0.0_dp),  'trace(pauli2) = 0')
  call assert(pauli3%trace() - (0.0_dp,0.0_dp),  'trace(pauli3) = 0')

  call subsection('Testing commutation relations:')
  call assert(pauli1*pauli2 - pauli2*pauli1 - (0.0_dp,2.0_dp)*pauli3, '[pauli1,pauli2] = 2i pauli3')
  call assert(pauli2*pauli3 - pauli3*pauli2 - (0.0_dp,2.0_dp)*pauli1, '[pauli2,pauli3] = 2i pauli1')
  call assert(pauli3*pauli1 - pauli1*pauli3 - (0.0_dp,2.0_dp)*pauli2, '[pauli3,pauli1] = 2i pauli2')



  ! Test the constructors and assignment operators in 'module_spin'
  call section('Construction and assignment of spin matrices')
  call subsection('Testing construction and assignment:')
  q%matrix(1,1) = (1.0_dp,2.0_dp)
  q%matrix(1,2) = (3.0_dp,4.0_dp)
  q%matrix(2,1) = (5.0_dp,6.0_dp)
  q%matrix(2,2) = (7.0_dp,8.0_dp)
  rv = [ 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp ]
  p  = spin(rv)
  call assert(p - q, 'Construction from a real vector')
  p  = rv
  call assert(p - q, 'Assignment from a real vector')
  rv = p
  call assert(rv - [ 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp ], 'Assignment to a real vector')
  cv = [ (1.0_dp,2.0_dp), (3.0_dp,4.0_dp), (5.0_dp,6.0_dp), (7.0_dp,8.0_dp) ]
  p  = spin(cv)
  call assert(p - q, 'Construction from a complex vector')
  p  = cv
  call assert(p - q, 'Assignment from a complex vector')
  cv = p
  call assert(cv - [ (1.0_dp, 2.0_dp), (3.0_dp, 4.0_dp), (5.0_dp, 6.0_dp), (7.0_dp, 8.0_dp) ], 'Assignment to a complex vector')
  p  = spin(q%matrix)
  call assert(p - q, 'Construction from a complex matrix')
  p  = q%matrix
  call assert(p - q, 'Assignment from a complex matrix')
  cm = q%matrix
  call assert([cm(1,1) - q%matrix(1,1), cm(1,2) - q%matrix(1,2), cm(2,1) - q%matrix(2,1), cm(2,2) - q%matrix(2,2)], &
              'Assignment to a complex matrix')
  p  = spin(q)
  call assert(p - q, 'Construction from a spin object')
  p  = q
  call assert(p - q, 'Assignment from a spin object')
  rs = 2.0_dp
  p  = spin(rs)
  call assert(p - 2.0_dp*pauli0, 'Construction from a real scalar')
  cs = (2.0_dp,1.0_dp)
  p  = spin(cs)
  call assert(p - (2.0_dp,1.0_dp)*pauli0, 'Construction from a complex scalar')



  ! Test the matrix algebra operators in 'module_spin'
  call section('Overloaded operators for matrix algebra')
  call subsection('Testing elementary matrix operations:')
  p = [ ( 1.00_dp, 1.00_dp), ( 1.00_dp,-1.00_dp), ( 2.00_dp, 3.00_dp), ( 2.00_dp,-3.00_dp) ]
  q = [ ( 1.00_dp, 2.00_dp), ( 2.00_dp,-1.00_dp), ( 3.00_dp,-3.00_dp), ( 2.00_dp, 2.00_dp) ]
  r = [ ( 2.00_dp, 3.00_dp), ( 3.00_dp,-2.00_dp), ( 5.00_dp, 0.00_dp), ( 4.00_dp,-1.00_dp) ]
  call assert(p+q - r, 'Matrix addition')
  r = [ ( 0.00_dp,-1.00_dp), (-1.00_dp, 0.00_dp), (-1.00_dp, 6.00_dp), ( 0.00_dp,-5.00_dp) ]
  call assert(p-q - r, 'Matrix subtraction')
  r = [ (-1.00_dp,-3.00_dp), ( 7.00_dp, 1.00_dp), (-7.00_dp,-8.00_dp), (17.00_dp, 2.00_dp) ]
  call assert(p*q - r, 'Matrix multiplication')
  r = [ ( 0.60_dp,-0.20_dp), ( 0.00_dp, 0.00_dp), ( 1.48_dp,-0.56_dp), (-0.20_dp, 0.00_dp) ]
  call assert((p .divr. q) - r, 'Matrix division (right)')
  r = [ (-3.50_dp,+4.00_dp), ( 4.00_dp,-1.50_dp), ( 3.50_dp,+5.00_dp), ( 0.00_dp,-3.50_dp) ]
  call assert((p .divl. q) - r, 'Matrix division (left)')
end program
