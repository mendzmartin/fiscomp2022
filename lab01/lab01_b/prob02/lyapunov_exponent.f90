program lyapunov_exponent
    use module_presition, only: pr => dp
    implicit none
    integer                                  :: i,j,N,Nt
    real(pr), allocatable                    :: xn(:),r(:)
    real(pr)                                 :: ra,rb,x0,delta,nn


    Nt = 300._pr
    N = 1500._pr

    ra = 2._pr
    rb = 4._pr

    x0 = 0.6
    nn = 1._pr/N
    r = ra + (rb-ra)*(/(i, i=0,1000)/)/1000._pr
    delta = 0._pr
    open(5,file='xn-r-lyapunov.d',status='replace')
    do i = 1,size(r)-1
        allocate(xn(N))
        xn(1) = x0
        do j = 1,size(xn)-1
            xn(j+1) = r(i)*xn(j)*(1-xn(j))
            if (j > Nt) then
                delta = delta + log(abs(r(i)-2._pr*r(i)*xn(j)))
            endif
        enddo
        write(5,*) r(i), delta*nn
        delta = 0._pr
        deallocate(xn)
    enddo
    close(5)



end program lyapunov_exponent