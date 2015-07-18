
      program test_materials
        use module_assert
        use module_conductor
        use module_superconductor
        integer              :: n
        real(dp)             :: erg(100) = [ ((n/50.0_dp), n=1,100) ]
        type(conductor)      :: a, b, c
        type(superconductor) :: s

        a = conductor(erg)
        b = conductor(erg)
        c = conductor(erg)
        s = superconductor(erg)
        call connect(a, b, 1/3.0_dp, 1/3.0_dp)
        call connect(c, s, 0.3_dp, 0.3_dp)
        call connect(s, c, 0.3_dp, 0.3_dp)

        do n=1,128
          print *,s%get_gap(s%location(n))
        end do
        call s%update
        
        call a%update

        open(unit=1, file='test_materials.dat', status='NEW')
        call a%write_dos(1)
        call b%write_dos(1)
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
