program test_spin
  use module_assert
  use module_spin

  call section('Calibration of assert functions')
  call assert(.true.,                'Logical true                          SUCCESS == ')
  call assert(.false.,               'Logical false                         FAILURE == ')
  call assert(0,                     'Integer 0                             SUCCESS == ')
  call assert(1,                     'Integer 1                             FAILURE == ')
  call assert(0.0_dp,                'Real 0                                SUCCESS == ')
  call assert(1.0_dp,                'Real 1                                FAILURE == ')
  call assert(spin(0.0_dp),          'Spin matrix 0                         SUCCESS == ')
  call assert(spin(1.0_dp),          'Spin matrix 1                         FAILURE == ')
  call assert(spin((1.0_dp,1.0_dp)), 'Spin matrix 1+i                       FAILURE == ')

  call section('Pauli matrices')
  call pauli0%print(' pauli0')
  call pauli1%print(' pauli1')
  call pauli2%print(' pauli2')
  call pauli3%print(' pauli3')

  call subsection('Testing norms:')
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
  

  !call assert(pauli0, pauli0, 'pauli0 == pauli0')
  !call assert(pauli0, pauli1, 'pauli0 == pauli1')


end program
