! Comparación entre el método implícito y la solución exacta
program heateq_comparison_02
    use module_precision
    use module_tridiag_matrix
    implicit none
    real(dp), allocatable   :: ux_old(:),ux_new(:)      ! vectores de temperatura (paso anterior y paso posterior)
    real(dp), allocatable   :: diag(:),diag_sup(:),diag_inf(:) ! diagonales de la matriz tridiagonal
    integer(sp), parameter  :: n=8_sp                   ! numero de elementos finitos dimension espacial
    real(dp), parameter     :: u0_Li=0._dp, u0_Lf=0._dp ! condiciones de borde
    real(dp)                :: u0_ti=100._dp            ! condicion inicial
    real(dp)                :: param_t                  ! tiempo caracteristico (paso temporal) [t]
    real(dp)                :: param_x                  ! longitud caracteristica (long total barra) [L]
    real(dp)                :: D                        ! thermal diffusivity (D=K/C*rho [L^2]/[t])
    real(dp)                :: alpha=1.E-7_dp           ! coeficiente discretización
    real(dp)                :: t_step_adim, x_step_adim ! deltatiempo y deltaespacio (adimensionales)
    integer(sp)             :: i,j,k,istat              ! variables de loop y variable de control
    real(dp), parameter     :: t_write_1_noadm=180._dp,t_write_2_noadm=1800._dp ! tiempos en dimensiones para escrituras
    integer(dp)             :: t_write_1,t_write_2      ! pasos temporales entre escrituras
    20 format (E8.2,x,E8.2)
    open( 10, file = '../results/result_02_aprox_implicit.dat', status = 'replace', action = 'write', iostat = istat )
    ! ANÁLISIS DIMENSIONAL
    D=237._dp*(1._dp/(900._dp*2700._dp))
    param_x = 1._dp
    param_t = param_x*param_x*(1._dp/D)
    x_step_adim=1._dp/(real(n,dp)+1._dp)
    t_step_adim=alpha*x_step_adim*x_step_adim
    ! CARGAMOS DIAGONALES CENTRAL, SUPERIOR E INFERIOR
    allocate(diag(n),diag_sup(n),diag_inf(n))
    diag_sup(n)=0._dp
    diag_inf(1)=0._dp
    do i=1,n
        diag(i)=1._dp+2._dp*alpha
        if (i/=1) diag_inf(i) = -alpha
        if (i/=n) diag_sup(i) = -alpha
    end do
    ! CARGAMOS DATOS INICIALES DE TEMPERATURA
    allocate(ux_old(n+2))
    ux_old(1)=u0_Li
    ux_old(2:n+1)=u0_ti
    ux_old(n+2)=u0_Lf
    ! ADIMENSIONALIZAMOS LOS TIEMPOS DE ESCRITURA
    t_write_1=int(t_write_1_noadm*(1._dp/(param_t*t_step_adim)),dp)
    t_write_2=int(t_write_2_noadm*(1._dp/(param_t*t_step_adim)),dp)
    ! APLICAMOS MÉTODO IMPLÍCITO
    allocate(ux_new(n+2))
    ux_new(1)=ux_old(1)
    ux_new(n+2)=ux_old(n+2)
    do j=1,(t_write_1-1)
        call implicit_method(n,diag,diag_sup,diag_inf,ux_old(2:n+1),ux_new(2:n+1))
        do k=2,(n+1); ux_old(k)=ux_new(k); end do
    end do
    ! ESCRIBIMOS VALORES LUEGO DE t_write_1 PASOS TEMPORALES
    call implicit_method(n,diag,diag_sup,diag_inf,ux_old(2:n+1),ux_new(2:n+1))
    do j=1,n+2; write(10,20) real(j-1,dp)*x_step_adim, ux_new(j); end do
    do j=t_write_1,(t_write_2-1)
        call implicit_method(n,diag,diag_sup,diag_inf,ux_old(2:n+1),ux_new(2:n+1))
        do k=2,(n+1); ux_old(k)=ux_new(k); end do
    end do
    ! ESCRIBIMOS VALORES LUEGO DE t_write_2 PASOS TEMPORALES
    call implicit_method(n,diag,diag_sup,diag_inf,ux_old(2:n+1),ux_new(2:n+1))
    do j=1,n+2; write(10,20) real(j-1,dp)*x_step_adim, ux_new(j); end do
    close(10)
    deallocate(diag,diag_sup,diag_inf,ux_old,ux_new)
end program heateq_comparison_02