! Problema 01
program random_generator
    use module_precision
    use module_random_generator
    use module_mzran
    use module_mt19937
    
    implicit none
    integer(sp), parameter   :: nrand=10000000_sp   ! cantidad de nros random
    real(dp),    parameter   :: lambda=1.5_dp       ! parametro de la distribución exp
    real(dp),    parameter   :: pi=4._dp*atan(1._dp)
    integer(sp)              :: seed                ! semilla
    integer(sp)              :: i                   ! loop variable
    real(dp),    allocatable :: x_rand0(:)          ! stocastic variable using ran0 random generator
    real(dp),    allocatable :: x_rand2(:)          ! stocastic variable using ran2 random generator
    real(dp),    allocatable :: x_mzran(:)          ! stocastic variable using mzran random generator
    real(dp),    allocatable :: x_mtwis(:)          ! stocastic variable using mt19937 random generator
    
    ! definimos semilla
    seed=121648_sp
    call sgrnd(seed) ! init mt19937 random generator

    allocate(x_rand0(nrand),x_rand2(nrand),x_mzran(nrand),x_mtwis(nrand))

    ! exponential distribution (x \in {0,Infinity})
    do i=1,nrand
        x_rand0(i)=-log(ran0(seed))*(1._dp/lambda)
        x_rand2(i)=-log(ran2(seed))*(1._dp/lambda)
        x_mzran(i)=-log(rmzran())*(1._dp/lambda)
        x_mtwis(i)=-log(real(grnd(),dp))*(1._dp/lambda)
    end do
    ! creamos histogramas
    call histogram('../results/histogram_exp_ran0.dat',x_rand0,nrand,1)
    call histogram('../results/histogram_exp_ran2.dat',x_rand2,nrand,1)
    call histogram('../results/histogram_exp_rmzran.dat',x_mzran,nrand,1)
    call histogram('../results/histogram_exp_mt19937.dat',x_mtwis,nrand,1)

    ! x^-2 distribution (x \in {1,Infinity})
    do i=1,nrand
        x_rand0(i)=1._dp/ran0(seed)
        x_rand2(i)=1._dp/ran2(seed)
        x_mzran(i)=1._dp/rmzran()
        x_mtwis(i)=1._dp/real(grnd(),dp)
    end do
    ! creamos histogramas
    call histogram('../results/histogram_pot_ran0.dat',x_rand0,nrand,2)
    call histogram('../results/histogram_pot_ran2.dat',x_rand2,nrand,2)
    call histogram('../results/histogram_pot_rmzran.dat',x_mzran,nrand,2)
    call histogram('../results/histogram_pot_mt19937.dat',x_mtwis,nrand,2)

    ! gaussian distribution (x \in {-Infinity,Infinity})
    do i=1,nrand
        x_rand0(i)=sqrt(-2*log(ran0(seed)))*cos(2*pi*ran0(seed))
        x_rand2(i)=sqrt(-2*log(ran2(seed)))*cos(2*pi*ran2(seed))
        x_mzran(i)=sqrt(-2*log(rmzran()))*cos(2*pi*rmzran())
        x_mtwis(i)=sqrt(-2*log(real(grnd(),dp)))*cos(2*pi*real(grnd(),dp))
    end do
    ! creamos histogramas
    call histogram('../results/histogram_gauss_ran0.dat',x_rand0,nrand,3)
    call histogram('../results/histogram_gauss_ran2.dat',x_rand2,nrand,3)
    call histogram('../results/histogram_gauss_rmzran.dat',x_mzran,nrand,3)
    call histogram('../results/histogram_gauss_mt19937.dat',x_mtwis,nrand,3)

    deallocate(x_rand0,x_rand2,x_mzran,x_mtwis)

end program random_generator

! subrutina para crear e imprmir histograma
! partiendo de un vector de datos normalizado
subroutine histogram(file_name,x_vector,x_dim,type_var)
    use module_precision

    implicit none
    character(len=*), intent(in) :: file_name
    integer(sp),      intent(in) :: x_dim           ! dimension
    real(dp),         intent(in) :: x_vector(x_dim) ! normalized data
    integer(sp),      intent(in) :: type_var        ! distribution type

    integer(sp), parameter   :: n_bins=100_sp ! numbers of bins (JUGAR CON ESTE VALOR)
    real(dp),    allocatable :: bins_ends(:)  ! bins ends vector
    integer(sp), allocatable :: counter(:)    ! counter vector of bins
    integer(sp)              :: max_value     ! maximun counter value
    integer(sp)              :: i,j,istat     ! loop and control variables

    ! armamos el vector de bins (normalizado)
    allocate(bins_ends(n_bins+1),counter(n_bins))

    select case(type_var)
    case(1) ! exponential distribution
        do i=1,n_bins+1; bins_ends(i)=real(i-1,dp)*(1._dp/real(n_bins-1,dp));end do
    case(2) ! pot distribution
        do i=1,n_bins+1; bins_ends(i)=1._dp+real(i-1,dp)*(1._dp/real(n_bins-1,dp));end do
    case(3) ! gauss distribution
        do i=1,n_bins+1; bins_ends(i)=-1._dp+2._dp*real(i-1,dp)*(1._dp/real(n_bins-1,dp));end do
    case default
        write(*,*) 'Invalid type variable'
    end select
    
    ! llenamos el vector contador de bins
    do i=1,n_bins
        counter(i)=0
        do j=1,x_dim
            if ((bins_ends(i)<=x_vector(j)).and.(bins_ends(i+1)>=x_vector(j))) then
                counter(i)=counter(i)+1
            endif
        enddo
    enddo

    max_value=maxval(counter(:))

    if (max_value==0) write(*,*) 'math error'

    ! escribimos el histograma
    open(10,file=file_name,status='replace',action='write',iostat=istat)
    20 format(E11.4,x,E11.4)
    write(*,*) 'istat=', istat
    do i = 1,n_bins; write(10,20) bins_ends(i),real(counter(i)*1._dp/max_value,dp); enddo
    close(10)
end subroutine histogram