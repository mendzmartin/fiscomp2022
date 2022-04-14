program lyapunov_exponent
    use module_presition
    implicit none
    integer(sp)            :: i,j,istat
    integer(sp), parameter :: N=1500_sp, Nt=300_sp, Nr=1000_sp
    real(dp), allocatable  :: xn(:)
    real(dp), parameter    :: r_min=2._dp,r_max=4._dp,x0 = 0.6
    real(dp)               :: r
    real(dp)               :: lambda_N

    open( 10, file = './lyapunov_exp.dat', status = 'replace', action = 'write', iostat = istat )
    write(*,*) istat
    allocate(xn(N))
    do i = 1,Nr
        r = r_min + abs(r_max-r_min)*(1._dp/(real(Nr,dp)-1._dp))*(real(i,dp)-1._dp)
        xn(1) = x0
        lambda_N = 0._dp
        do j = 2,N
            xn(j) = r*xn(j-1_sp)*(1._dp-xn(j-1_sp))
            if (j > Nt) then
                lambda_N = lambda_N + log(abs(r*(1._dp-2._dp*xn(j-1_sp))))
            endif
        enddo
        write(10,*) r, lambda_N*(1._dp/real(N,dp))
    enddo
    deallocate(xn)
    close(10)

end program lyapunov_exponent