
      program test_materials
        use module_assert
        use module_conductor
        use module_superconductor
        integer              :: n
        real(dp)             :: erg(100) = [ ((n/50.0_dp), n=1,100) ] !, (2.0_dp + n/(cosh(5.0_dp)-2.0_dp), n=1,100) ]
        type(conductor)      :: m
        type(superconductor) :: s

        m = conductor(erg)
        s = superconductor(erg)
        call connect(m, s, 0.3_dp, 0.3_dp)

        do n=1,128
          print *,s%get_gap(s%location(n))
        end do
        call m%update
        call s%update
        do n=1,128
          print *,s%get_gap(s%location(n))
        end do

        open(unit=1, file='test_materials.dat', status='NEW')
        call m%write_dos(1)
        call s%write_dos(1)
        close(unit=1)


        !call b%update
        !print *,b%state(50,64)%get_dos()
        !call c%update
        !print *,c%state(50,64)%get_dos()
        !call b%update
        !print *,b%state(50,64)%get_dos()
        !call a%update
        !print *,a%state(50,64)%get_dos()
        !call b%update
        !print *,b%state(50,64)%get_dos()

        !print *, metal%state(50,64)%get_dos()

        !system = [ metal, metal, metal ]

        !subroutine write_test(fd)
        !end subroutine
      end program
