program orbits_diagram_chaos
    use module_presition
    implicit none
    integer(sp)             :: i,j,N,istat
    real(dp), allocatable   :: xn(:)
    real(dp)                :: r, r_min,r_max
    integer(sp), parameter  :: Nr=1000_sp, n_trans=300_sp
    real(dp), parameter     :: x0=0.6_dp

    r_min = 3.4_dp
    r_max = 4._dp

    open( 10, file = './bif_graph_01.dat', status = 'replace', action = 'write', iostat = istat )
    write(*,*) istat
    N = 600_sp
    allocate(xn(N))
    do i = 1,Nr
        r = r_min + abs(r_max-r_min)*(1._dp/(real(Nr,dp)-1._dp))*(real(i,dp)-1._dp)
        xn(1) = x0
        do j = 2,N
            xn(j) = r*xn(j-1_sp)*(1._dp-xn(j-1_sp))  ! ecuacion logistica
            if (j > n_trans) then
                write(10,*) r, xn(j)
            endif
        enddo
    enddo
    deallocate(xn)
    close(10)

    r_max = 3.847_dp
    r_min = 3.8568_dp

    open( 11, file = './bif_graph_02.dat', status = 'replace', action = 'write', iostat = istat )
    write(*,*) istat
    N = 1500_sp
    allocate(xn(N))
    do i = 1,Nr
        r = r_min + abs(r_max-r_min)*(1._dp/(real(Nr,dp)-1._dp))*(real(i,dp)-1._dp)
        xn(1) = x0
        do j = 2,N
            xn(j) = r*xn(j-1_sp)*(1._dp-xn(j-1_sp))  ! ecuacion logistica
            if (j > n_trans) then
                write(11,*) r, xn(j)
            endif
        enddo
    enddo
    deallocate(xn)
    close(11)

end program orbits_diagram_chaos