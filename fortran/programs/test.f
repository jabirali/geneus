! This program is used to test the fundamental data structures and procedures.
!
! Author:  Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
! Created: 2015-07-13
! Updated: 2015-08-09

program test_core
  use mod_assert
  use mod_spin
  use mod_green

  real(wp)    :: rs, rv(8), rw(32)
  complex(wp) :: cs, cv(4), cm(2,2)
  type(spin)  :: p, q, r
  type(green) :: s0, s1, s2
  integer     :: n


  ! Calibrate the test procedures in 'module_assert'
  call section('Calibration of assert subroutines')
  call subsection('Subroutines operating on a single value:')
  call assert(.true.,                'Logical true                                     SUCCESS == ')
  call assert(.false.,               'Logical false                                    FAILURE == ')
  call assert(0,                     'Integer 0                                        SUCCESS == ')
  call assert(1,                     'Integer 1                                        FAILURE == ')
  call assert(0.0_wp,                'Real    0                                        SUCCESS == ')
  call assert(1.0_wp,                'Real    1                                        FAILURE == ')
  call assert((0.0_wp,0.0_wp),       'Complex 0                                        SUCCESS == ')
  call assert((0.0_wp,1.0_wp),       'Complex i                                        FAILURE == ')
  call assert(spin(0.0_wp),          'Spin matrix 0                                    SUCCESS == ')
  call assert(spin(1.0_wp),          'Spin matrix 1                                    FAILURE == ')
  call assert(spin((0.0_wp,1.0_wp)), 'Spin matrix i                                    FAILURE == ')

  call subsection('Subroutines operating on arrays:')
  call assert([.true., .true.],                  'Logical [true, true ]                            SUCCESS == ')
  call assert([.true., .false.],                 'Logical [true, false]                            FAILURE == ')
  call assert([.false.,.false.],                 'Logical [false,false]                            FAILURE == ')
  call assert([0, 0],                            'Integer [0,0]                                    SUCCESS == ')
  call assert([1, 0],                            'Integer [1,0]                                    FAILURE == ')
  call assert([1, 1],                            'Integer [1,1]                                    FAILURE == ')
  call assert([0.0_wp,0.0_wp],                   'Real    [0,0]                                    SUCCESS == ')
  call assert([1.0_wp,0.0_wp],                   'Real    [1,0]                                    FAILURE == ')
  call assert([1.0_wp,1.0_wp],                   'Real    [1,1]                                    FAILURE == ')
  call assert([(0.0_wp,0.0_wp),(0.0_wp,0.0_wp)], 'Complex [0,0]                                    SUCCESS == ')
  call assert([(0.0_wp,1.0_wp),(0.0_wp,0.0_wp)], 'Complex [i,0]                                    FAILURE == ')
  call assert([(0.0_wp,1.0_wp),(0.0_wp,1.0_wp)], 'Complex [i,i]                                    FAILURE == ')



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
  call assert(pauli0%trace() - (2.0_wp,0.0_wp),  'trace(pauli0) = 2')
  call assert(pauli1%trace() - (0.0_wp,0.0_wp),  'trace(pauli1) = 0')
  call assert(pauli2%trace() - (0.0_wp,0.0_wp),  'trace(pauli2) = 0')
  call assert(pauli3%trace() - (0.0_wp,0.0_wp),  'trace(pauli3) = 0')

  call subsection('Testing commutation relations:')
  call assert(pauli1*pauli2 - pauli2*pauli1 - (0.0_wp,2.0_wp)*pauli3, '[pauli1,pauli2] = 2i pauli3')
  call assert(pauli2*pauli3 - pauli3*pauli2 - (0.0_wp,2.0_wp)*pauli1, '[pauli2,pauli3] = 2i pauli1')
  call assert(pauli3*pauli1 - pauli1*pauli3 - (0.0_wp,2.0_wp)*pauli2, '[pauli3,pauli1] = 2i pauli2')



  ! Test the constructors and assignment operators in 'module_spin'
  call section('Spin: Construction and assignment')
  call subsection('Testing construction and assignment:')
  q%matrix(1,1) = (1.0_wp,2.0_wp)
  q%matrix(1,2) = (3.0_wp,4.0_wp)
  q%matrix(2,1) = (5.0_wp,6.0_wp)
  q%matrix(2,2) = (7.0_wp,8.0_wp)
  rv = [ 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp ]
  p  = spin(rv)
  call assert(p - q, 'Construction from a real vector')
  p  = rv
  call assert(p - q, 'Assignment from a real vector')
  rv = p
  call assert(rv - [ 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp ], 'Assignment to a real vector')
  cv = [ (1.0_wp,2.0_wp), (3.0_wp,4.0_wp), (5.0_wp,6.0_wp), (7.0_wp,8.0_wp) ]
  p  = spin(cv)
  call assert(p - q, 'Construction from a complex vector')
  p  = cv
  call assert(p - q, 'Assignment from a complex vector')
  cv = p
  call assert(cv - [ (1.0_wp, 2.0_wp), (3.0_wp, 4.0_wp), (5.0_wp, 6.0_wp), (7.0_wp, 8.0_wp) ], 'Assignment to a complex vector')
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
  rs = 2.0_wp
  p  = spin(rs)
  call assert(p - 2.0_wp*pauli0, 'Construction from a real scalar')
  cs = (2.0_wp,1.0_wp)
  p  = spin(cs)
  call assert(p - (2.0_wp,1.0_wp)*pauli0, 'Construction from a complex scalar')



  ! Test the matrix algebra operators in 'module_spin'
  call section('Spin: Arithmetic operators and methods')
  call subsection('Testing matrix operations with scalars:')
  p = [ ( 1.00_wp, 1.00_wp), ( 1.00_wp,-1.00_wp), ( 2.00_wp, 3.00_wp), ( 2.00_wp,- 3.00_wp) ]
  q = [ ( 1.00_wp, 2.00_wp), ( 2.00_wp,-1.00_wp), ( 3.00_wp,-3.00_wp), ( 2.00_wp,  2.00_wp) ]
  r = [ ( 5.00_wp, 3.00_wp), ( 1.00_wp,-5.00_wp), (12.00_wp, 5.00_wp), ( 0.00_wp,-11.00_wp) ]
  call assert(p**2  - r, 'Exponentiation by an integer')
  r = [ ( 2.00_wp, 2.00_wp), ( 2.00_wp,-2.00_wp), ( 4.00_wp, 6.00_wp), ( 4.00_wp,-6.00_wp) ]
  call assert(2.0_wp*p - r, 'Multiplication by a real scalar (left)')
  call assert(p*2.0_wp - r, 'Multiplication by a real scalar (right)')
  r = [ (-2.00_wp, 2.00_wp), ( 2.00_wp, 2.00_wp), (-6.00_wp, 4.00_wp), ( 6.00_wp, 4.00_wp) ]
  call assert((0.0_wp,2.0_wp)*p - r, 'Multiplication by a complex scalar (left)')
  call assert(p*(0.0_wp,2.0_wp) - r, 'Multiplication by a complex scalar (right)')
  r = [ ( 0.50_wp, 0.50_wp), ( 0.50_wp,-0.50_wp), ( 1.00_wp, 1.50_wp), ( 1.00_wp,-1.50_wp) ]
  call assert(p/2.0_wp - r, 'Division by a real scalar')
  r = [ ( 0.50_wp,-0.50_wp), (-0.50_wp,-0.50_wp), ( 1.50_wp,-1.00_wp), (-1.50_wp,-1.00_wp) ]
  call assert(p/(0.0_wp,2.0_wp) - r, 'Division by a complex scalar')

  call subsection('Testing matrix operations with matrices:')
  r = [ ( 2.00_wp, 3.00_wp), ( 3.00_wp,-2.00_wp), ( 5.00_wp, 0.00_wp), ( 4.00_wp,-1.00_wp) ]
  call assert(p+q - r, 'Matrix addition')
  r = [ ( 0.00_wp,-1.00_wp), (-1.00_wp, 0.00_wp), (-1.00_wp, 6.00_wp), ( 0.00_wp,-5.00_wp) ]
  call assert(p-q - r, 'Matrix subtraction')
  r = [ (-1.00_wp,-3.00_wp), ( 7.00_wp, 1.00_wp), (-7.00_wp,-8.00_wp), (17.00_wp, 2.00_wp) ]
  call assert(p*q - r, 'Matrix multiplication')
  r = [ ( 0.60_wp,-0.20_wp), ( 0.00_wp, 0.00_wp), ( 1.48_wp,-0.56_wp), (-0.20_wp, 0.00_wp) ]
  call assert((p .divr. q) - r, 'Matrix division (right)')
  r = [ (-3.50_wp, 4.00_wp), ( 4.00_wp,-1.50_wp), ( 3.50_wp, 5.00_wp), ( 0.00_wp,-3.50_wp) ]
  call assert((p .divl. q) - r, 'Matrix division (left)')

  call subsection('Testing matrix methods:')
  r = [ ( 1.0000_wp, 2.0000_wp), ( 3.0000_wp, 4.0000_wp), ( 5.0000_wp, 6.0000_wp), ( 7.0000_wp, 8.0000_wp) ]
  q = [ (-0.5000_wp, 0.4375_wp), ( 0.2500_wp,-0.1875_wp), ( 0.3750_wp,-0.3125_wp), (-0.1250_wp, 0.0625_wp) ]
  call assert(r%inv()   - q,                     'Inverse:')
  call assert(r%trace() - (8.0_wp,10.0_wp),      'Trace:')
  call assert(r%norm()  - 14.282856857085701_wp, 'Norm:')
  call assert(r%min()   -  2.236067977499790_wp, 'Minimum:')
  call assert(r%max()   - 10.630145812734650_wp, 'Maximum:')


  ! Test construction and assignment in 'module_spin'
  call section('State: Construction and assignment')
  call subsection('Testing output:')
  s0 = green(pauli0, pauli1, pauli2, pauli3)
  call s0%print

  call subsection('Testing assignment operator:')
  rw = [ (n, n=1,32) ]
  s0 = rw
  s1 = green( spin([( 1.000000000000000_wp, 2.000000000000000_wp), ( 3.000000000000000_wp, 4.000000000000000_wp),   &
                    ( 5.000000000000000_wp, 6.000000000000000_wp), ( 7.000000000000000_wp, 8.000000000000000_wp)]), &
              spin([( 9.000000000000000_wp,10.000000000000000_wp), (11.000000000000000_wp,12.000000000000000_wp),   &
                    (13.000000000000000_wp,14.000000000000000_wp), (15.000000000000000_wp,16.000000000000000_wp)]), &
              spin([(17.000000000000000_wp,18.000000000000000_wp), (19.000000000000000_wp,20.000000000000000_wp),   &
                    (21.000000000000000_wp,22.000000000000000_wp), (23.000000000000000_wp,24.000000000000000_wp)]), &
              spin([(25.000000000000000_wp,26.000000000000000_wp), (27.000000000000000_wp,28.000000000000000_wp),   &
                    (29.000000000000000_wp,30.000000000000000_wp), (31.000000000000000_wp,32.000000000000000_wp)]))
  call assert(s0%g         - s1%g,                     'Importing from a real vector (gamma )')
  call assert(s0%gt        - s1%gt,                    'Importing from a real vector (gamma~)')
  call assert(s0%dg        - s1%dg,                    'Importing from a real vector (dgamma / dz)')
  call assert(s0%dgt       - s1%dgt,                   'Importing from a real vector (dgamma~/ dz)')
  rw = s0
  call assert(rw - [ (n, n=1,32) ],                    'Exporting to a real vector')

  call subsection('Testing construction of a BCS state with E=0 (gap real):')
  s0 = green( (0.0_wp,0.001_wp), (1.0_wp,0.0_wp) )
  s1 = green( spin([( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp,-0.999000499999875_wp),   &
                    ( 0.000000000000000_wp, 0.999000499999875_wp), ( 0.000000000000000_wp, 0.000000000000000_wp)]), &
              spin([( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp, 0.999000499999875_wp),   &
                    ( 0.000000000000000_wp,-0.999000499999875_wp), ( 0.000000000000000_wp, 0.000000000000000_wp)]), &
              spin([( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp, 0.000000000000000_wp),   &
                    ( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp, 0.000000000000000_wp)]), &
              spin([( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp, 0.000000000000000_wp),   &
                    ( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp, 0.000000000000000_wp)]))
  call assert(s0%g         - s1%g,                     'Riccati parameter gamma ')
  call assert(s0%gt        - s1%gt,                    'Riccati parameter gamma~')
  call assert(s0%dg        - s1%dg,                    'Derivative dgamma / dz')
  call assert(s0%dgt       - s1%dgt,                   'Derivative dgamma~/ dz')
  call assert(s0%get_dos() - 9.999995000002913e-04_wp, 'Density of states')

  call subsection('Testing construction of a BCS state with E=1 (gap real):')
  s0 = green( (1.0_wp,0.001_wp), (1.0_wp,0.0_wp) )
  s1 = green( spin([( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.968385128104009_wp,-0.030630683283800_wp),   &
                    (-0.968385128104009_wp, 0.030630683283800_wp), ( 0.000000000000000_wp, 0.000000000000000_wp)]), &
              spin([( 0.000000000000000_wp, 0.000000000000000_wp), (-0.968385128104009_wp, 0.030630683283800_wp),   &
                    ( 0.968385128104009_wp,-0.030630683283800_wp), ( 0.000000000000000_wp, 0.000000000000000_wp)]), &
              spin([( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp, 0.000000000000000_wp),   &
                    ( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp, 0.000000000000000_wp)]), &
              spin([( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp, 0.000000000000000_wp),   &
                    ( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp, 0.000000000000000_wp)]))
  call assert(s0%g         - s1%g,                     'Riccati parameter gamma ')
  call assert(s0%gt        - s1%gt,                    'Riccati parameter gamma~')
  call assert(s0%dg        - s1%dg,                    'Derivative dgamma / dz')
  call assert(s0%dgt       - s1%dgt,                   'Derivative dgamma~/ dz')
  call assert(s0%get_dos() - 15.823249311731509_wp,    'Density of states')

  call subsection('Testing construction of a BCS state with E=1 (gap imag):')
  s0 = green( (1.0_wp,0.001_wp), (0.0_wp,1.0_wp) )
  s1 = green( spin([( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.030630683283800_wp, 0.968385128104009_wp),   &
                    (-0.030630683283800_wp,-0.968385128104009_wp), ( 0.000000000000000_wp, 0.000000000000000_wp)]), &
              spin([( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.030630683283800_wp, 0.968385128104009_wp),   &
                    (-0.030630683283800_wp,-0.968385128104009_wp), ( 0.000000000000000_wp, 0.000000000000000_wp)]), &
              spin([( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp, 0.000000000000000_wp),   &
                    ( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp, 0.000000000000000_wp)]), &
              spin([( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp, 0.000000000000000_wp),   &
                    ( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp, 0.000000000000000_wp)]))
  call assert(s0%g         - s1%g,                     'Riccati parameter gamma ')
  call assert(s0%gt        - s1%gt,                    'Riccati parameter gamma~')
  call assert(s0%dg        - s1%dg,                    'Derivative dgamma / dz')
  call assert(s0%dgt       - s1%dgt,                   'Derivative dgamma~/ dz')
  call assert(s0%get_dos() - 15.823249311731509_wp,    'Density of states')

  call subsection('Testing construction of a BCS state with E=2 (gap real):')
  s0 = green( (2.0_wp,0.001_wp), (1.0_wp,0.0_wp) )
  s1 = green( spin([( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.267949096206123_wp,-0.000154700474229_wp),   &
                    (-0.267949096206123_wp, 0.000154700474229_wp), ( 0.000000000000000_wp, 0.000000000000000_wp)]), &
              spin([( 0.000000000000000_wp, 0.000000000000000_wp), (-0.267949096206123_wp, 0.000154700474229_wp),   &
                    ( 0.267949096206123_wp,-0.000154700474229_wp), ( 0.000000000000000_wp, 0.000000000000000_wp)]), &
              spin([( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp, 0.000000000000000_wp),   &
                    ( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp, 0.000000000000000_wp)]), &
              spin([( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp, 0.000000000000000_wp),   &
                    ( 0.000000000000000_wp, 0.000000000000000_wp), ( 0.000000000000000_wp, 0.000000000000000_wp)]))
  call assert(s0%g         - s1%g,                     'Riccati parameter gamma ')
  call assert(s0%gt        - s1%gt,                    'Riccati parameter gamma~')
  call assert(s0%dg        - s1%dg,                    'Derivative dgamma / dz')
  call assert(s0%dgt       - s1%dgt,                   'Derivative dgamma~/ dz')
  call assert(s0%get_dos() - 1.154700345929331_wp,     'Density of states')

  call section("State: Calculation of Green's functions")
  s0%g  = (0.05_wp,0.10_wp)*pauli0 + (0.15_wp,0.20_wp)*pauli1 + (0.25_wp,0.30_wp)*pauli2 + (0.35_wp,0.40_wp)*pauli3
  s0%gt = (0.07_wp,0.08_wp)*pauli0 + (0.17_wp,0.18_wp)*pauli1 + (0.27_wp,0.28_wp)*pauli2 + (0.37_wp,0.38_wp)*pauli3
  p = s0%get_g()
  q = [ ( 0.405054515669019_wp, 0.780842410470802_wp), (-0.007173683883726_wp, 0.148657413387066_wp),&
        (-0.128531295454375_wp,-0.092454113122009_wp), ( 0.697757267837770_wp, 0.639249425699371_wp) ]
  call assert(p - q, 'Extraction of g :')
  p = s0%get_gt()
  q = [ ( 0.396085295417136_wp, 0.788328900154937_wp), (-0.031115883503879_wp, 0.138205462567435_wp),&
        (-0.122527536337988_wp,-0.067029182934108_wp), ( 0.706726488089653_wp, 0.631762936015236_wp) ]
  call assert(p - q, 'Extraction of g~:')
  p = s0%get_f()
  q = [ ( 0.105780817590586_wp, 0.989337452267094_wp), ( 0.718065981755837_wp, 0.238681240077408_wp),&
        (-0.547511293361128_wp, 0.566856063696100_wp), (-0.380014141252089_wp,-0.736279794193328_wp) ]
  call assert(p - q, 'Extraction of f :')
  p = s0%get_ft()
  q = [ ( 0.192876524942350_wp, 0.959859203500891_wp), ( 0.749984387773050_wp, 0.196973454763113_wp),&
        (-0.495111926735743_wp, 0.612677489472185_wp), (-0.383659083438489_wp,-0.720682481281396_wp) ]
  call assert(p - q, 'Extraction of f~:')
  call assert(s0%get_f_s()  - (0.632788637558482_wp,-0.164087411809346_wp), 'Extraction of the singlet component of f :')
  call assert(s0%get_ft_s() - (0.622548157254396_wp,-0.207852017354536_wp), 'Extraction of the singlet component of f~:')
  call assert(s0%get_f_t()  - [(-0.242897479421337_wp,-0.862808623230211_wp), &
                               ( 0.126528829036883_wp, 0.137116661830751_wp), &
                               ( 0.085277344197354_wp, 0.402768651886754_wp)],&
              'Extraction of the triplet component of f :')
  call assert(s0%get_ft_t() - [(-0.288267804190420_wp,-0.840270842391144_wp), &
                               ( 0.119588361109748_wp, 0.095391279248070_wp), &
                               ( 0.127436230518653_wp, 0.404825472117649_wp)],&
              'Extraction of the triplet component of f~:')
  call assert(s0%get_f_ts([1.0_wp, 2.0_wp, 3.0_wp])  - [( 0.018999443660321_wp, 0.044266475435111_wp), &
                                                        ( 0.037998887320642_wp, 0.088532950870222_wp), &
                                                        ( 0.056998330980962_wp, 0.132799426305333_wp)],&
              'Extraction of the short-range triplet component of f :')
  call assert(s0%get_ft_ts([1.0_wp, 2.0_wp, 3.0_wp]) - [( 0.023801257827502_wp, 0.040356295175567_wp), &
                                                        ( 0.047602515655005_wp, 0.080712590351135_wp), &
                                                        ( 0.071403773482507_wp, 0.121068885526702_wp)],&
              'Extraction of the short-range triplet component of f~:')
  call assert(s0%get_f_tl([1.0_wp, 2.0_wp, 3.0_wp])  - [(-0.261896923081658_wp,-0.907075098665322_wp), &
                                                        ( 0.088529941716241_wp, 0.048583710960530_wp), &
                                                        ( 0.028279013216392_wp, 0.269969225581421_wp)],&
              'Extraction of the long-range triplet component of f :')
  call assert(s0%get_ft_tl([1.0_wp, 2.0_wp, 3.0_wp]) - [(-0.312069062017922_wp,-0.880627137566711_wp), &
                                                        ( 0.071985845454743_wp, 0.014678688896935_wp), &
                                                        ( 0.056032457036146_wp, 0.283756586590947_wp)],&
              'Extraction of the long-range triplet component of f~:')
  call assert(s0%get_dos() - 0.551405891753394_wp, 'Extraction of the density of states:')
contains
 subroutine section(msg)
   character(*),  intent(in) :: msg
   character(68)             :: str

   str = trim(msg)

   print *,''
   print *,'╒═══════════════════════════════════',&
           '═══════════════════════════════════╕'
   print *,'│  ',          str,               '│'
   print *,'╘═══════════════════════════════════',&
           '═══════════════════════════════════╛'
 end subroutine

 subroutine subsection(msg)
   character(*), intent(in) :: msg

   print *,' '
   print *,msg
 end subroutine
end program
