program orbits_diagram_chaos
    use module_presition
    implicit none
    integer(sp)             :: i,j,N,n_trans, istat
    real(dp), allocatable   :: xn(:),r(:)
    real(dp)                :: r_min,r_max,x0
    integer(sp), parameter  :: Nr=1000_sp

    r_min = 3.4_dp
    r_max = 4._dp

    x0 = 0.6_dp
    allocate(r(Nr))
    r = r_min + abs(r_max-r_min)*real((/(i, i=0,Nr)/),dp)*(1._dp/(real(Nr,dp)))

    open( 10, file = './bif_graph_01.dat', status = 'replace', action = 'write', iostat = istat )
    write(*,*) istat
    N = 600_sp
    n_trans = 300_sp
    do i = 1,(size(r)-1)
        allocate(xn(N))
        xn(1) = x0
        do j = 1,(size(xn)-1_sp)
            xn(j+1) = r(i)*xn(j)*(1._dp-xn(j))  ! ecuacion logistica
            if (j > n_trans) then
                write(10,*) r(i), xn(j)
            endif
        enddo
        deallocate(xn)
    enddo
    close(10)

    r_max = 3.847_dp
    r_min = 3.8568_dp

    r = r_min + abs(r_max-r_min)*real((/(i, i=0,Nr)/),dp)*(1._dp/(real(Nr,dp)))

    open( 10, file = './bif_graph_02.dat', status = 'replace', action = 'write', iostat = istat )
    write(*,*) istat
    N = 1500_sp
    do i = 1,(size(r)-1)
        allocate(xn(N))
        xn(1) = x0
        do j = 1,(size(xn)-1_sp)
            xn(j+1) = r(i)*xn(j)*(1-xn(j))
            if (j > n_trans) then
                write(10,*) r(i), xn(j)
            endif
        enddo
        deallocate(xn)
    enddo
    close(10)

    deallocate(r)

end program orbits_diagram_chaos