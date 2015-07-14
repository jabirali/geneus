! This program is used to test the data structures and procedures declared in 'modules/*.f90'.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-13
! Updated: 2015-07-13

program test_core
  use module_assert
  use module_spin
  use module_state

  real(dp)    :: rs, rv(8), rw(32)
  complex(dp) :: cs, cv(4), cm(2,2)
  type(spin)  :: p, q, r
  type(state) :: s0, s1, s2


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
  call section('Spin: Pauli matrices')
  call subsection('Testing output:')
  call pauli0%print('pauli0')
  call pauli1%print('pauli1')
  call pauli2%print('pauli2')
  call pauli3%print('pauli3')

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
  call section('Spin: Construction and assignment')
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
  call section('Spin: Arithmetic operators and methods')
  call subsection('Testing matrix operations with scalars:')
  p = [ ( 1.00_dp, 1.00_dp), ( 1.00_dp,-1.00_dp), ( 2.00_dp, 3.00_dp), ( 2.00_dp,- 3.00_dp) ]
  q = [ ( 1.00_dp, 2.00_dp), ( 2.00_dp,-1.00_dp), ( 3.00_dp,-3.00_dp), ( 2.00_dp,  2.00_dp) ]
  r = [ ( 5.00_dp, 3.00_dp), ( 1.00_dp,-5.00_dp), (12.00_dp, 5.00_dp), ( 0.00_dp,-11.00_dp) ]
  call assert(p**2  - r, 'Exponentiation by an integer')
  r = [ ( 2.00_dp, 2.00_dp), ( 2.00_dp,-2.00_dp), ( 4.00_dp, 6.00_dp), ( 4.00_dp,-6.00_dp) ]
  call assert(2.0_dp*p - r, 'Multiplication by a real scalar (left)')
  call assert(p*2.0_dp - r, 'Multiplication by a real scalar (right)')
  r = [ (-2.00_dp, 2.00_dp), ( 2.00_dp, 2.00_dp), (-6.00_dp, 4.00_dp), ( 6.00_dp, 4.00_dp) ]
  call assert((0.0_dp,2.0_dp)*p - r, 'Multiplication by a complex scalar (left)')
  call assert(p*(0.0_dp,2.0_dp) - r, 'Multiplication by a complex scalar (right)')
  r = [ ( 0.50_dp, 0.50_dp), ( 0.50_dp,-0.50_dp), ( 1.00_dp, 1.50_dp), ( 1.00_dp,-1.50_dp) ]
  call assert(p/2.0_dp - r, 'Division by a real scalar')
  r = [ ( 0.50_dp,-0.50_dp), (-0.50_dp,-0.50_dp), ( 1.50_dp,-1.00_dp), (-1.50_dp,-1.00_dp) ]
  call assert(p/(0.0_dp,2.0_dp) - r, 'Division by a complex scalar')

  call subsection('Testing matrix operations with matrices:')
  r = [ ( 2.00_dp, 3.00_dp), ( 3.00_dp,-2.00_dp), ( 5.00_dp, 0.00_dp), ( 4.00_dp,-1.00_dp) ]
  call assert(p+q - r, 'Matrix addition')
  r = [ ( 0.00_dp,-1.00_dp), (-1.00_dp, 0.00_dp), (-1.00_dp, 6.00_dp), ( 0.00_dp,-5.00_dp) ]
  call assert(p-q - r, 'Matrix subtraction')
  r = [ (-1.00_dp,-3.00_dp), ( 7.00_dp, 1.00_dp), (-7.00_dp,-8.00_dp), (17.00_dp, 2.00_dp) ]
  call assert(p*q - r, 'Matrix multiplication')
  r = [ ( 0.60_dp,-0.20_dp), ( 0.00_dp, 0.00_dp), ( 1.48_dp,-0.56_dp), (-0.20_dp, 0.00_dp) ]
  call assert((p .divr. q) - r, 'Matrix division (right)')
  r = [ (-3.50_dp, 4.00_dp), ( 4.00_dp,-1.50_dp), ( 3.50_dp, 5.00_dp), ( 0.00_dp,-3.50_dp) ]
  call assert((p .divl. q) - r, 'Matrix division (left)')

  call subsection('Testing matrix methods:')
  r = [ ( 1.0000_dp, 2.0000_dp), ( 3.0000_dp, 4.0000_dp), ( 5.0000_dp, 6.0000_dp), ( 7.0000_dp, 8.0000_dp) ]
  q = [ (-0.5000_dp, 0.4375_dp), ( 0.2500_dp,-0.1875_dp), ( 0.3750_dp,-0.3125_dp), (-0.1250_dp, 0.0625_dp) ]
  call assert(r%inv()   - q,                     'Inverse:')
  call assert(r%trace() - (8.0_dp,10.0_dp),      'Trace:')
  call assert(r%norm()  - 14.282856857085701_dp, 'Norm:')
  call assert(r%min()   -  2.236067977499790_dp, 'Minimum:')
  call assert(r%max()   - 10.630145812734650_dp, 'Maximum:')


  ! Test construction and assignment in 'module_spin'
  call section('State: Construction and assignment')
  call subsection('Testing output:')
  s0 = state(pauli0, pauli1, pauli2, pauli3)
  call s0%print

  call subsection('Testing assignment operator:')
  rw = [ (n, n=1,32) ]
  s0 = rw
  s1 = state( spin([( 1.000000000000000_dp, 2.000000000000000_dp), ( 3.000000000000000_dp, 4.000000000000000_dp),   &
                    ( 5.000000000000000_dp, 6.000000000000000_dp), ( 7.000000000000000_dp, 8.000000000000000_dp)]), &
              spin([( 9.000000000000000_dp,10.000000000000000_dp), (11.000000000000000_dp,12.000000000000000_dp),   &
                    (13.000000000000000_dp,14.000000000000000_dp), (15.000000000000000_dp,16.000000000000000_dp)]), &
              spin([(17.000000000000000_dp,18.000000000000000_dp), (19.000000000000000_dp,20.000000000000000_dp),   &
                    (21.000000000000000_dp,22.000000000000000_dp), (23.000000000000000_dp,24.000000000000000_dp)]), &
              spin([(25.000000000000000_dp,26.000000000000000_dp), (27.000000000000000_dp,28.000000000000000_dp),   &
                    (29.000000000000000_dp,30.000000000000000_dp), (31.000000000000000_dp,32.000000000000000_dp)]))
  call assert(s0%g         - s1%g,                     'Importing from a real vector (γ )')
  call assert(s0%gt        - s1%gt,                    'Importing from a real vector (γ~)')
  call assert(s0%dg        - s1%dg,                    'Importing from a real vector (dγ / dz)')
  call assert(s0%dgt       - s1%dgt,                   'Importing from a real vector (dγ~/ dz)')
  rw = s0
  call assert(rw - [ (n, n=1,32) ],                    'Exporting to a real vector　')

  call subsection('Testing construction of a BCS state with ϵ=0 and Δ=1:')
  s0 = state( (0.0_dp,0.001_dp), (1.0_dp,0.0_dp) )
  s1 = state( spin([( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp,-0.999000499999875_dp),   &
                    ( 0.000000000000000_dp, 0.999000499999875_dp), ( 0.000000000000000_dp, 0.000000000000000_dp)]), &
              spin([( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp, 0.999000499999875_dp),   &
                    ( 0.000000000000000_dp,-0.999000499999875_dp), ( 0.000000000000000_dp, 0.000000000000000_dp)]), &
              spin([( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp, 0.000000000000000_dp),   &
                    ( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp, 0.000000000000000_dp)]), &
              spin([( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp, 0.000000000000000_dp),   &
                    ( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp, 0.000000000000000_dp)]))
  call assert(s0%g         - s1%g,                     'Riccati parameter γ ')
  call assert(s0%gt        - s1%gt,                    'Riccati parameter γ~')
  call assert(s0%dg        - s1%dg,                    'Derivative dγ / dz')
  call assert(s0%dgt       - s1%dgt,                   'Derivative dγ~/ dz')
  call assert(s0%get_dos() - 9.999995000002913e-04_dp, 'Density of states　')

  call subsection('Testing construction of a BCS state with ϵ=1 and Δ=1:')
  s0 = state( (1.0_dp,0.001_dp), (1.0_dp,0.0_dp) )
  s1 = state( spin([( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.968385128104009_dp,-0.030630683283800_dp),   &
                    (-0.968385128104009_dp, 0.030630683283800_dp), ( 0.000000000000000_dp, 0.000000000000000_dp)]), &
              spin([( 0.000000000000000_dp, 0.000000000000000_dp), (-0.968385128104009_dp, 0.030630683283800_dp),   &
                    ( 0.968385128104009_dp,-0.030630683283800_dp), ( 0.000000000000000_dp, 0.000000000000000_dp)]), &
              spin([( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp, 0.000000000000000_dp),   &
                    ( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp, 0.000000000000000_dp)]), &
              spin([( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp, 0.000000000000000_dp),   &
                    ( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp, 0.000000000000000_dp)]))
  call assert(s0%g         - s1%g,                     'Riccati parameter γ ')
  call assert(s0%gt        - s1%gt,                    'Riccati parameter γ~')
  call assert(s0%dg        - s1%dg,                    'Derivative dγ / dz')
  call assert(s0%dgt       - s1%dgt,                   'Derivative dγ~/ dz')
  call assert(s0%get_dos() - 15.823249311731509_dp,    'Density of states　')

  call subsection('Testing construction of a BCS state with ϵ=1 and Δ=i:')
  s0 = state( (1.0_dp,0.001_dp), (0.0_dp,1.0_dp) )
  s1 = state( spin([( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.030630683283800_dp, 0.968385128104009_dp),   &
                    (-0.030630683283800_dp,-0.968385128104009_dp), ( 0.000000000000000_dp, 0.000000000000000_dp)]), &
              spin([( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.030630683283800_dp, 0.968385128104009_dp),   &
                    (-0.030630683283800_dp,-0.968385128104009_dp), ( 0.000000000000000_dp, 0.000000000000000_dp)]), &
              spin([( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp, 0.000000000000000_dp),   &
                    ( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp, 0.000000000000000_dp)]), &
              spin([( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp, 0.000000000000000_dp),   &
                    ( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp, 0.000000000000000_dp)]))
  call assert(s0%g         - s1%g,                     'Riccati parameter γ ')
  call assert(s0%gt        - s1%gt,                    'Riccati parameter γ~')
  call assert(s0%dg        - s1%dg,                    'Derivative dγ / dz')
  call assert(s0%dgt       - s1%dgt,                   'Derivative dγ~/ dz')
  call assert(s0%get_dos() - 15.823249311731509_dp,    'Density of states　')

  call subsection('Testing construction of a BCS state with ϵ=2 and Δ=1:')
  s0 = state( (2.0_dp,0.001_dp), (1.0_dp,0.0_dp) )
  s1 = state( spin([( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.267949096206123_dp,-0.000154700474229_dp),   &
                    (-0.267949096206123_dp, 0.000154700474229_dp), ( 0.000000000000000_dp, 0.000000000000000_dp)]), &
              spin([( 0.000000000000000_dp, 0.000000000000000_dp), (-0.267949096206123_dp, 0.000154700474229_dp),   &
                    ( 0.267949096206123_dp,-0.000154700474229_dp), ( 0.000000000000000_dp, 0.000000000000000_dp)]), &
              spin([( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp, 0.000000000000000_dp),   &
                    ( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp, 0.000000000000000_dp)]), &
              spin([( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp, 0.000000000000000_dp),   &
                    ( 0.000000000000000_dp, 0.000000000000000_dp), ( 0.000000000000000_dp, 0.000000000000000_dp)]))
  call assert(s0%g         - s1%g,                     'Riccati parameter γ ')
  call assert(s0%gt        - s1%gt,                    'Riccati parameter γ~')
  call assert(s0%dg        - s1%dg,                    'Derivative dγ / dz')
  call assert(s0%dgt       - s1%dgt,                   'Derivative dγ~/ dz')
  call assert(s0%get_dos() - 1.154700345929331_dp,     'Density of states　')

  call section("State: Calculation of Green's functions")
  s0%g  = (0.05_dp,0.10_dp)*pauli0 + (0.15_dp,0.20_dp)*pauli1 + (0.25_dp,0.30_dp)*pauli2 + (0.35_dp,0.40_dp)*pauli3
  s0%gt = (0.07_dp,0.08_dp)*pauli0 + (0.17_dp,0.18_dp)*pauli1 + (0.27_dp,0.28_dp)*pauli2 + (0.37_dp,0.38_dp)*pauli3
  p = s0%get_g()
  q = [ ( 0.405054515669019_dp, 0.780842410470802_dp), (-0.007173683883726_dp, 0.148657413387066_dp),&
        (-0.128531295454375_dp,-0.092454113122009_dp), ( 0.697757267837770_dp, 0.639249425699371_dp) ]
  call assert(p - q, 'Extraction of g :')
  p = s0%get_gt()
  q = [ ( 0.396085295417136_dp, 0.788328900154937_dp), (-0.031115883503879_dp, 0.138205462567435_dp),&
        (-0.122527536337988_dp,-0.067029182934108_dp), ( 0.706726488089653_dp, 0.631762936015236_dp) ]
  call assert(p - q, 'Extraction of g~:')
  p = s0%get_f()
  q = [ ( 0.105780817590586_dp, 0.989337452267094_dp), ( 0.718065981755837_dp, 0.238681240077408_dp),&
        (-0.547511293361128_dp, 0.566856063696100_dp), (-0.380014141252089_dp,-0.736279794193328_dp) ]
  call assert(p - q, 'Extraction of f :')
  p = s0%get_ft()
  q = [ ( 0.192876524942350_dp, 0.959859203500891_dp), ( 0.749984387773050_dp, 0.196973454763113_dp),&
        (-0.495111926735743_dp, 0.612677489472185_dp), (-0.383659083438489_dp,-0.720682481281396_dp) ]
  call assert(p - q, 'Extraction of f~:')
  call assert(s0%get_f_s()  - (0.632788637558482_dp,-0.164087411809346_dp), 'Extraction of the singlet component of f :')
  call assert(s0%get_ft_s() - (0.622548157254396_dp,-0.207852017354536_dp), 'Extraction of the singlet component of f~:')
  call assert(s0%get_f_t()  - [(-0.242897479421337_dp,-0.862808623230211_dp), &
                               ( 0.126528829036883_dp, 0.137116661830751_dp), &
                               ( 0.085277344197354_dp, 0.402768651886754_dp)],&
              'Extraction of the triplet component of f :')
  call assert(s0%get_ft_t() - [(-0.288267804190420_dp,-0.840270842391144_dp), &
                               ( 0.119588361109748_dp, 0.095391279248070_dp), &
                               ( 0.127436230518653_dp, 0.404825472117649_dp)],&
              'Extraction of the triplet component of f~:')
  call assert(s0%get_f_ts([1.0_dp, 2.0_dp, 3.0_dp])  - [( 0.018999443660321_dp, 0.044266475435111_dp), &
                                                        ( 0.037998887320642_dp, 0.088532950870222_dp), &
                                                        ( 0.056998330980962_dp, 0.132799426305333_dp)],&
              'Extraction of the short-range triplet component of f :')
  call assert(s0%get_ft_ts([1.0_dp, 2.0_dp, 3.0_dp]) - [( 0.023801257827502_dp, 0.040356295175567_dp), &
                                                        ( 0.047602515655005_dp, 0.080712590351135_dp), &
                                                        ( 0.071403773482507_dp, 0.121068885526702_dp)],&
              'Extraction of the short-range triplet component of f~:')
  call assert(s0%get_f_tl([1.0_dp, 2.0_dp, 3.0_dp])  - [(-0.261896923081658_dp,-0.907075098665322_dp), &
                                                        ( 0.088529941716241_dp, 0.048583710960530_dp), &
                                                        ( 0.028279013216392_dp, 0.269969225581421_dp)],&
              'Extraction of the long-range triplet component of f :')
  call assert(s0%get_ft_tl([1.0_dp, 2.0_dp, 3.0_dp]) - [(-0.312069062017922_dp,-0.880627137566711_dp), &
                                                        ( 0.071985845454743_dp, 0.014678688896935_dp), &
                                                        ( 0.056032457036146_dp, 0.283756586590947_dp)],&
              'Extraction of the long-range triplet component of f~:')
  call assert(s0%get_dos() - 0.551405891753394_dp, 'Extraction of the density of states:')
end program
