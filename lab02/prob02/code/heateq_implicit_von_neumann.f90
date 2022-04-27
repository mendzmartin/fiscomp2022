program heateq_implicit_von_neumann
    use module_precision;use module_tridiag_matrix
    implicit none
    real(dp),    parameter   :: u0_Li=0._dp, u0_Lf=0._dp                 ! condiciones de borde
    real(dp),    parameter   :: t_step_adim=0.001_dp,x_step_adim=0.05_dp ! deltaT y deltaX adimensionales
    integer(sp), parameter   :: t_write=2_sp                             ! pasos temporales entre escrituras
    real(dp)                 :: pi=4._dp*atan(1.0)
    real(dp)                 :: param_t,param_x
    real(dp)                 :: D                        ! thermal diffusivity
    real(dp)                 :: alpha
    integer(sp)              :: n,a,b                    ! numero de elementos finitos y variables auxiliares
    integer(sp)              :: i,j,istat
    real(dp),    allocatable :: ux_old(:),ux_new(:),diag(:),diag_sup(:),diag_inf(:)
    20 format(E13.6,x,E13.6)
    open(10,file='../results/result_01_implicit_vn.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat
    ! parametros para adimensionalizar
    D=237._dp*(1._dp/(900._dp*2700._dp)) ! D = K/C*rho [L^2]/[t]
    param_x=1._dp                        ! longitud caracteristica (long total barra en metros) [L]
    param_t=param_x*param_x*(1._dp/D)    ! tiempo caracteristico (paso temporal en segundos) [t]
    alpha=t_step_adim*(1._dp/(x_step_adim*x_step_adim))
    ! calculamos la parte entera del modulo y del resto del inverso de x_step_adim
    a=int(mod(1._dp,x_step_adim),sp);b=int(1._dp/x_step_adim);n=(a+b-1_sp)*(1_sp/(1_sp-a))
    ! cargamos diagonales central, superior e inferior
    allocate(diag(n+2),diag_sup(n+2),diag_inf(n+2))
    diag_sup(n+2)=0._dp;diag_inf(1)=0._dp
    do i=1,n+2
        diag(i)=1._dp+2._dp*alpha
        if (i/=1) diag_inf(i)=-alpha;if (i/=n+2) diag_sup(i)=-alpha
    end do
    ! cambiamos los extremos de las diagonales iferior y superior
    diag_inf(n+1)=diag_inf(n+1)*2._dp;diag_sup(2)=diag_sup(2)*2._dp
    ! cargamos datos iniciales de temperatura
    allocate(ux_old(n+2))
    do i=1,n+2
        ux_old(i)=cos(pi*real(i-1,dp)*x_step_adim)
        write(10,20) real(i-1,dp)*x_step_adim,ux_old(i)
    end do
    ! aplicamos método implicito
    allocate(ux_new(n+2))
    do i=1,10
        do j=1,(t_write-1)
            ux_old(1)=ux_old(1)+2._dp*alpha*x_step_adim*u0_Li
            ux_old(n+2)=ux_old(n+2)-2._dp*alpha*x_step_adim*u0_Lf
            call implicit_method(n+2,diag,diag_sup,diag_inf,ux_old(:),ux_new(:))
            ux_old(:)=ux_new(:)
        end do
        call implicit_method(n+2,diag,diag_sup,diag_inf,ux_old(:),ux_new(:))
        do j=1,n+2;write(10,20) real(j-1,dp)*x_step_adim,ux_new(j);end do
    end do
    close(10);deallocate(diag,diag_sup,diag_inf,ux_old,ux_new)
end program heateq_implicit_von_neumann