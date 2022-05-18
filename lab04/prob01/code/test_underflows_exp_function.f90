! make clean && make test_underflows_exp_function.o && ./test_underflows_exp_function.o
program test_underflows_exp_function
    use module_precision
    implicit none
    real(dp)    :: x,result
    integer(sp) :: i
    ! tiny(x) returns the smallest positive (non zero) number in the model of the type of x
    write(*,*) 'abs(log(tiny(x)))=',abs(log(tiny(1._dp)))
    do i=1,800
        x=real(i,dp)
        if (x<abs(log(tiny(result)))) then; result=exp(-x)
        else; result=0._dp; end if
        write(*,*) x,result
    end do
end program test_underflows_exp_function
! https://gcc.gnu.org/onlinedocs/gcc-4.5.4/gfortran/TINY.html