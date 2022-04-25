program heateq_comparison_03
    use module_precision
    use module_tridiag_matrix
    implicit none
    real(dp), allocatable   :: ux_old(:),ux_new(:)      ! vectores de temperatura (paso anterior y paso posterior)
    real(dp), allocatable   :: diag(:),diag_sup(:),diag_inf(:) ! diagonales de la matriz tridiagonal
    real(dp),    allocatable :: B_matrix(:,:),A_matrix(:,:)
    integer(sp), parameter  :: n=8_sp                   ! numero de elementos finitos dimension espacial
    real(dp), parameter     :: u0_Li=0._dp, u0_Lf=0._dp ! condiciones de borde
    real(dp)                :: u0_ti=100._dp            ! condicion inicial
    real(dp)                :: param_t                  ! tiempo caracteristico (paso temporal) [t]
    real(dp)                :: param_x                  ! longitud caracteristica (long total barra) [L]
    real(dp)                :: D                        ! thermal diffusivity (D=K/C*rho [L^2]/[t])
    real(dp)                :: alpha=1.E-7_dp           ! coeficiente discretización
    real(dp)                :: t_step_adim, x_step_adim ! deltatiempo y deltaespacio (adimensionales)
    integer(sp)             :: i,istat              ! variables de loop y variable de control
    real(dp), parameter     :: t_write_1_noadm=180._dp,t_write_2_noadm=1800._dp ! tiempos en dimensiones para escrituras
    integer(dp)             :: t_write_1,t_write_2      ! pasos temporales entre escrituras
    20 format (E8.2,x,E8.2)
    open( 10, file = '../results/result_02_aprox_cranknicolson.dat', status = 'replace', action = 'write', iostat = istat )
    ! ANÁLISIS DIMENSIONAL
    D=237._dp*(1._dp/(900._dp*2700._dp))
    param_x = 1._dp
    param_t = param_x*param_x*(1._dp/D)
    x_step_adim=1._dp/(real(n,dp)+1._dp)
    t_step_adim=alpha*x_step_adim*x_step_adim
    ! ADIMENSIONALIZAMOS LOS TIEMPOS DE ESCRITURA
    t_write_1=int(t_write_1_noadm*(1._dp/(param_t*t_step_adim)),dp)
    t_write_2=int(t_write_2_noadm*(1._dp/(param_t*t_step_adim)),dp)
    ! CARGAMOS DIAGONALES CENTRAL, SUPERIOR E INFERIOR
    allocate(diag(n),diag_sup(n),diag_inf(n))
    diag_sup(n)=0._dp
    diag_inf(1)=0._dp
    do i=1,n
        diag(i)=2._dp*(1._dp/alpha+1._dp)
        if (i/=1) diag_inf(i)=-1._dp
        if (i/=n) diag_sup(i)=-1._dp
    end do
    ! CARGAMOS MATRIZ CUADRADA PARA USAR MÉTODO IMPLÍCITO
    allocate(B_matrix(n,n))
    B_matrix=0._dp ! elementos nulos fuera de la tribanda
    do i=1,n
        B_matrix(i,i)=2._dp*(1._dp/alpha-1._dp) ! diagonal principal
        if (i/=1) B_matrix(i,i-1)=1._dp         ! diagonal inferior
        if (i/=n) B_matrix(i,i+1)=1._dp         ! diagonal superior
    end do
    ! CARGAMOS DATOS INICIALES DE TEMPERATURA
    allocate(A_matrix(n,1)) ! matriz auxiliar p/matmul (column vector rank=2)
    allocate(ux_old(n+2))
    ux_old(1)=u0_Li
    A_matrix(1:n,1)=u0_ti
    ux_old(n+2)=u0_Lf
    !CARGAMOS DATOS INICIALES PARA APLCIAR MÉTODO IMPLÍCITO
    A_matrix=matmul(B_matrix,A_matrix)
    ux_old(2:n+1)=A_matrix(1:n,1)
    ! APLICAMOS MÉTODO IMPLÍCITO
    allocate(ux_new(n+2))
    ux_new(1)=ux_old(1)
    ux_new(n+2)=ux_old(n+2)
    do i=1,(t_write_1-1)
        call implicit_method(n,diag,diag_sup,diag_inf,data_vector=ux_old(2:n+1),unknown_vector=ux_new(2:n+1))
        A_matrix(1:n,1)=ux_new(2:n+1)
        A_matrix=matmul(B_matrix,A_matrix)
        ux_old(2:n+1)=A_matrix(1:n,1)
    end do
    ! ESCRIBIMOS VALORES LUEGO DE t_write_1 PASOS TEMPORALES
    call implicit_method(n,diag,diag_sup,diag_inf,ux_old(2:n+1),ux_new(2:n+1))
    do i=1,n+2
        write(10,20) real(i-1,dp)*x_step_adim, ux_new(i)
    end do
    do i=t_write_1,(t_write_2-1)
        call implicit_method(n,diag,diag_sup,diag_inf,data_vector=ux_old(2:n+1),unknown_vector=ux_new(2:n+1))
        A_matrix(1:n,1)=ux_new(2:n+1)
        A_matrix=matmul(B_matrix,A_matrix)
        ux_old(2:n+1)=A_matrix(1:n,1)
    end do
    ! ESCRIBIMOS VALORES LUEGO DE t_write_2 PASOS TEMPORALES
    call implicit_method(n,diag,diag_sup,diag_inf,ux_old(2:n+1),ux_new(2:n+1))
    do i=1,n+2
        write(10,20) real(i-1,dp)*x_step_adim, ux_new(i)
    end do
    close(10)
    deallocate(diag,diag_sup,diag_inf,ux_old,ux_new)
    deallocate(B_matrix,A_matrix)
end program heateq_comparison_03