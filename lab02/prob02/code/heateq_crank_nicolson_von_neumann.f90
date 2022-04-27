program heateq_crank_nicolson_von_neumann
    use module_precision;use module_tridiag_matrix
    implicit none
    real(dp),    parameter   :: u0_Li=0._dp, u0_Lf=0._dp                 ! condiciones de borde de Von Neumann (t=t0)
    real(dp),    parameter   :: u1_Li=0._dp, u1_Lf=0._dp                 ! condiciones de borde de Von Neumann (t=t1)
    real(dp),    parameter   :: t_step_adim=0.001_dp,x_step_adim=0.05_dp ! deltaT y deltaX adimensionales
    integer(sp), parameter   :: t_write=2_sp                             ! pasos temporales entre escrituras
    real(dp),    parameter   :: pi=4._dp*atan(1.0)
    real(dp),    parameter   :: D=237._dp*(1._dp/(900._dp*2700._dp)) ! thermal diffusivity(D=(K/C*rho)[L^2]/[t])
    real(dp),    parameter   :: param_x=1._dp                     ! longitud caracteristica (long total barra en metros) [L]
    real(dp),    parameter   :: param_t=param_x*param_x*(1._dp/D) ! tiempo caracteristico (paso temporal en segundos) [t]
    real(dp),    parameter   :: alpha=t_step_adim*(1._dp/(x_step_adim*x_step_adim))
    integer(sp)              :: n,a,b                    ! numero de elementos finitos y variables auxiliares
    integer(sp)              :: i,j,istat
    real(dp),    allocatable :: ux_old(:),ux_new(:),diag(:),diag_sup(:),diag_inf(:)
    real(dp),    allocatable :: B_matrix(:,:),A_matrix(:,:)
    20 format(E13.6,x,E13.6,x,E13.6)
    open(10,file='../results/result_01_cranknicolson_vn.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat
    write(*,*) 't_adim = ',10*t_step_adim
    ! calculamos la parte entera del modulo y del resto del inverso de x_step_adim
    a=int(mod(1._dp,x_step_adim),sp);b=int(1._dp/x_step_adim);n=(a+b-1_sp)*(1_sp/(1_sp-a))
    ! CARGAMOS DIAGONALES CENTRAL, SUPERIOR E INFERIOR
    allocate(diag(n+2),diag_sup(n+2),diag_inf(n+2))
    diag(:)=2._dp*(1._dp/alpha+1._dp);diag_inf(:)=-1._dp;diag_sup(:)=-1._dp
    ! cambiamos los extremos de las diagonales iferior y superior
    diag_inf(n+2)=diag_inf(n+2)*2._dp;diag_sup(1)=diag_sup(1)*2._dp
    diag_sup(n+2)=0._dp;diag_inf(1)=0._dp
    ! CARGAMOS MATRIZ CUADRADA PARA USAR MÉTODO IMPLÍCITO
    allocate(B_matrix(n+2,n+2))
    B_matrix=0._dp ! elementos nulos fuera de la tribanda
    do i=1,n+2
        B_matrix(i,i)=2._dp*(1._dp/alpha-1._dp) ! diagonal principal
        if (i/=1) B_matrix(i,i-1)=1._dp         ! diagonal inferior
        if (i/=n+2) B_matrix(i,i+1)=1._dp       ! diagonal superior
    end do
    ! cambiamos los extremos de las diagonales iferior y superior de la matriz B
    B_matrix(n+2,n+1)=B_matrix(n+2,n+1)*2._dp;B_matrix(1,2)=B_matrix(1,2)*2._dp
    ! CARGAMOS DATOS INICIALES DE TEMPERATURA
    allocate(A_matrix(n+2,1)) ! matriz auxiliar p/matmul (column vector rank=2)
    allocate(ux_old(n+2))
    do i=1,n+2
        ux_old(i)=cos(pi*real(i-1,dp)*x_step_adim);A_matrix(i,1)=ux_old(i)
        write(10,20) 0.0_dp,real(i-1,dp)*x_step_adim,ux_old(i)
    end do
    write(10,*) ''
    ! CARGAMOS DATOS INICIALES PARA APLCIAR MÉTODO IMPLÍCITO
    A_matrix=matmul(B_matrix,A_matrix);ux_old(:)=A_matrix(:,1)
    ! APLICAMOS MÉTODO IMPLÍCITO
    allocate(ux_new(n+2))
    do i=1,1000
        do j=1,(t_write-1)
            ux_old(1)=ux_old(1)+2._dp*alpha*x_step_adim*(u0_Li-u1_Lf)
            ux_old(n+2)=ux_old(n+2)-2._dp*alpha*x_step_adim*(u0_Lf-u1_Lf)
            call implicit_method(n+2,diag,diag_sup,diag_inf,data_vector=ux_old(:),unknown_vector=ux_new(:))
            A_matrix(:,1)=ux_new(:);A_matrix=matmul(B_matrix,A_matrix);ux_old(:)=A_matrix(:,1)
        end do
        call implicit_method(n+2,diag,diag_sup,diag_inf,ux_old(:),ux_new(:))
        do j=1,n+2;write(10,20) t_write*i*t_step_adim*param_t,real(j-1,dp)*x_step_adim, ux_new(j);end do
        write(10,*) ''
    end do
    close(10);deallocate(diag,diag_sup,diag_inf,ux_old,ux_new,B_matrix,A_matrix)
end program heateq_crank_nicolson_von_neumann