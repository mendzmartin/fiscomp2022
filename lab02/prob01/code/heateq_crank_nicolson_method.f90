program heateq_crank_nicolson_method
    use module_precision
    use module_tridiag_matrix
    implicit none
    real(dp)                 :: param_t,param_x
    real(dp)                 :: D                        ! thermal diffusivity
    real(dp)                 :: alpha
    real(dp)                 :: t_step_adim,x_step_adim
    integer(sp)              :: i,j,istat
    real(dp),    allocatable :: ux_old(:),ux_new(:),diag(:),diag_sup(:),diag_inf(:)
    real(dp),    allocatable :: B_matrix(:,:),A_matrix(:,:)
    integer(sp), parameter   :: n=98_sp                  ! numero de elementos finitos
    real(dp),    parameter   :: u0_Li=0._dp, u0_Lf=0._dp ! condiciones de borde
    real(dp),    parameter   :: u0_ti=100._dp            ! condicion inicial
    integer(sp), parameter   :: t_write=300_sp           ! pasos temporales entre escrituras
    20 format(E8.2,x,E8.2)
    open(10,file='../results/result_01_cranknicolson.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat
    ! ANALISIS DIMENSIONAL
    D=237._dp*(1._dp/(900._dp*2700._dp))    ! D = K/C*rho [L^2]/[t]
    param_x=1._dp                           ! longitud caracteristica (long total barra) [L]
    param_t=param_x*param_x*(1._dp/D)       ! tiempo caracteristico (paso temporal) [t]
    t_step_adim=0.3_dp*(1._dp/param_t)
    x_step_adim=1._dp/(real(n,dp)+1._dp)
    write(*,*) 'x_step_adim = ',x_step_adim
    alpha=t_step_adim*(1._dp/(x_step_adim*x_step_adim))
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
    write(10,20) 0._dp, ux_old(1)
    do i=2,n+1
        A_matrix(i-1,1)=u0_ti
        write(10,20) real(i-1,dp)*x_step_adim, A_matrix(i-1,1)
    end do
    ux_old(n+2)=u0_Lf
    write(10,20) real(n+1,dp)*x_step_adim, ux_old(n+2)
    ! CARGAMOS DATOS INICIALES PARA APLCIAR MÉTODO IMPLÍCITO
    A_matrix=matmul(B_matrix,A_matrix)
    ux_old(2:n+1)=A_matrix(1:n,1)
    ! APLICAMOS MÉTODO IMPLÍCITO
    allocate(ux_new(n+2))
    ux_new(1)=ux_old(1)
    ux_new(n+2)=ux_old(n+2)
    do i=1,50
        do j=1,(t_write-1)
            call implicit_method(n,diag,diag_sup,diag_inf,data_vector=ux_old(2:n+1),unknown_vector=ux_new(2:n+1))
            A_matrix(1:n,1)=ux_new(2:n+1)
            A_matrix=matmul(B_matrix,A_matrix)
            ux_old(2:n+1)=A_matrix(1:n,1)
        end do
        call implicit_method(n,diag,diag_sup,diag_inf,ux_old(2:n+1),ux_new(2:n+1))
        do j=1,n+2
            write(10,20) real(j-1,dp)*x_step_adim, ux_new(j)
        end do
    end do
    close(10)
    ! CONTROLAMOS TEMALIZACIÓN EN EL CENTRO DE LA BARRA
    write(*,'(A20,E10.4)') 'T(L/2,t_final) = ', ux_new(50)
    deallocate(diag,diag_sup,diag_inf,ux_old,ux_new)
    deallocate(B_matrix,A_matrix)
end program heateq_crank_nicolson_method