type, public :: spinactive
  real(wp)                 :: conductance   = 0.0      !! Interfacial conductance
  real(wp)                 :: polarization  = 0.0      !! Interfacial spin-polarization
  real(wp)                 :: spinmixing    = 0.0      !! Interfacial 1st-order spin-mixing
  real(wp)                 :: secondorder   = 0.0      !! Interfacial 2nd-order spin-mixing
  real(wp), dimension(1:3) :: magnetization = [0,0,1]  !! Interfacial magnetization direction
  real(wp), dimension(1:3) :: misalignment  = [0,0,0]  !! Interfacial magnetization misalignment

  complex(wp), dimension(1:4,1:4) :: M  = 0.0          !! Magnetization matrix (transmission)
  complex(wp), dimension(1:4,1:4) :: M0 = 0.0          !! Magnetization matrix (reflection, this  side)
  complex(wp), dimension(1:4,1:4) :: M1 = 0.0          !! Magnetization matrix (reflection, other side)
end type
