! make clean && make test_pointers.o && ./test_pointers.o

program test_pointers
    use module_precision
    implicit none

    real(dp), target :: r,s
    real(dp), pointer :: p1,p2

    r=5.3_dp; s=8.1_dp
    p1 => r ! guarda en p1 la dirección de r y otra información descriptiva
    p2 => p1
    print *,p2+4._dp,p1,r
    p1 => s
    print *,p2+4._dp,p1,r
    ! p2+4._dp=9.3_dp, p1=5.3_dp,r=5.3_dp
    ! p2+4._dp=9.3_dp, p1=8.1_dp,r=5.3_dp
end program test_pointers