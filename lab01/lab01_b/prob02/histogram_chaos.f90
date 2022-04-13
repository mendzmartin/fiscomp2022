program histogram_chaos
    use module_presition

    implicit none

    integer(sp), parameter  :: n_trans=300  ! discarded steps number in transitory regime
    integer(sp), parameter  :: n_stat=3 ! used steps number in steady state regime
    integer(sp), parameter  :: n_bins=100   ! numbers of bins (intervals number)
    integer(sp)             :: i,j,istat
    real(dp), parameter     :: r=4._dp,x0=0.6_dp ! parameter and initial condition
    real(dp)                :: x_min,x_max ! min and max values of logistical equation
    real(dp), allocatable   :: x_stat(:) ! logistical vector in steady state regime
    real(dp), allocatable   :: x(:) ! logitical vector
    real(dp), allocatable   :: bins(:) ! bins vector
    real(sp), allocatable   :: counter(:)   ! counter vector of bins

    allocate(x(n_trans+n_stat))

    ! llenamos el vector (ecuación logística)
    x(1) = x0
    do i = 2_sp,size(x)
        x(i) = r*x(i-1)*(1._dp-x(i-1))
    enddo

    ! borramos el regimen transitorio
    allocate(x_stat(n_stat))
    do i = 1,n_stat
        if (i>n_trans) then
            x_stat(i-n_trans) = x(i)
        endif
    enddo
    deallocate(x)

    ! normalizamos la ecuación logistica
    x_min = minval(x_stat)
    x_max = maxval(x_stat)
    x_stat = (x_stat-x_min)*(1._dp/(x_max-x_min))

    ! armamos el vector de bins
    allocate(bins(n_bins),counter(n_bins))
    bins = real((/(i, i=0,n_bins)/),dp)*(1._dp/real(n_bins,dp))

    ! llenamos el vector contador de bins
    do i = 1_sp,(size(bins)-1_sp)
        counter(i) = 0_sp
        do j = 1,size(x_stat)
            if ( (bins(i) <= x_stat(j)) .and. (bins(i+1) >= x_stat(j)) ) then
                counter(i) = counter(i) + 1_sp
            endif
        enddo
    enddo

    deallocate(x_stat)

    open( 10, file = './histogram.dat', status = 'replace', action = 'write', iostat = istat )
    write(*,*) istat
    do i = 1,(size(bins)-1_sp)
        write(10,*) bins(i), counter(i)
    enddo

    deallocate(bins, counter)

end program histogram_chaos