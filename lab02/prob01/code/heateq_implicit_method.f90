program heateq_implicit_method
    use module_precision
    use module_tridiag_matrix
    implicit none
    integer(sp), parameter   :: n=98_sp                  ! numero de elementos finitos
    real(dp),    parameter   :: u0_Li=0._dp, u0_Lf=0._dp ! condiciones de borde
    real(dp),    parameter   :: u0_ti=100._dp            ! condicion inicial
    integer(sp), parameter   :: t_write=300_sp           ! pasos temporales entre escrituras
    real(dp),    allocatable :: ux_old(:),ux_new(:),diag(:),diag_sup(:),diag_inf(:)
    real(dp)                 :: param_t,param_x
    real(dp)                 :: D                        ! thermal diffusivity
    real(dp)                 :: alpha
    real(dp)                 :: t_step_adim,x_step_adim
    integer(sp)              :: i,j,k,istat
    20 format(E8.2,x,E8.2)
    open(10,file='../results/result_01_implicit.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat
    ! parametros para adimensionalizar
    D=237._dp*(1._dp/(900._dp*2700._dp))        ! D = K/C*rho [L^2]/[t]
    param_x=1._dp                              ! longitud caracteristica (long total barra) [L]
    param_t=param_x*param_x*(1._dp/D)          ! tiempo caracteristico (paso temporal) [t]
    t_step_adim=0.3_dp*(1._dp/param_t)
    x_step_adim=1._dp/(real(n,dp)+1._dp)
    alpha=t_step_adim*(1._dp/(x_step_adim*x_step_adim))
    ! cargamos diagonales central, superior e inferior
    allocate(diag(n),diag_sup(n),diag_inf(n))
    diag_sup(n)=0._dp
    diag_inf(1)=0._dp
    do i=1,n
        diag(i)=1._dp+2._dp*alpha
        if (i/=1) diag_inf(i)=-alpha
        if (i/=n) diag_sup(i)=-alpha
    end do
    ! cargamos datos iniciales de temperatura
    allocate(ux_old(n+2))
    ux_old(1)=u0_Li
    write(10,20) 0._dp,ux_old(1)
    do i=2,n+1
        ux_old(i)=u0_ti
        write(10,20) real(i-1,dp)*x_step_adim,ux_old(i)
    end do
    ux_old(n+2)=u0_Lf
    write(10,20) real(n+1,dp)*x_step_adim,ux_old(n+2)
    ! aplicamos método implicito
    allocate(ux_new(n+2))
    ux_new(1)=ux_old(1)
    ux_new(n+2)=ux_old(n+2)
    do i=1,50
        do j=1,(t_write-1)
            call implicit_method(n,diag,diag_sup,diag_inf,ux_old(2:n+1),ux_new(2:n+1))
            ux_old(2:n+1)=ux_new(2:n+1)
        end do
        call implicit_method(n,diag,diag_sup,diag_inf,ux_old(2:n+1),ux_new(2:n+1))
        do j=1,n+2
            write(10,20) real(j-1,dp)*x_step_adim,ux_new(j)
        end do
    end do
    close(10)
    deallocate(diag,diag_sup,diag_inf,ux_old,ux_new)
end program heateq_implicit_method