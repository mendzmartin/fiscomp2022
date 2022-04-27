program heateq_comparison_02
    use module_precision;use module_tridiag_matrix
    implicit none
    real(dp),    parameter   :: u0_Li=0._dp, u0_Lf=0._dp                            ! condiciones de borde (von neumann)
    real(dp),    parameter   :: t_step_adim=0.001_dp,x_step_adim=0.05_dp            ! deltaT y deltaX adimensionales
    real(dp),    parameter   :: t_write_adim=1._dp                                  ! tiempo adimensional de evolución
    real(dp),    parameter    :: pi=4._dp*atan(1.0)
    real(dp),    parameter   :: D=237._dp*(1._dp/(900._dp*2700._dp))                ! thermal diffusivity(D=(K/C*rho)[L^2]/[t])
    real(dp),    parameter   :: param_x=1._dp                                       ! longitud caracteristica (long total barra en metros) [L]
    real(dp),    parameter   :: param_t=param_x*param_x*(1._dp/D)                   ! tiempo caracteristico (paso temporal en segundos) [t]
    real(dp),    parameter   :: alpha=t_step_adim*(1._dp/(x_step_adim*x_step_adim))
    integer(sp), parameter   :: t_write=int(t_write_adim*(1._dp/t_step_adim),sp)    ! pasos temporales de evolución
    integer(sp)              :: n,a,b                                               ! numero de elementos finitos y variables auxiliares
    integer(sp)              :: i,j,istat
    real(dp)                 :: T_exact                                             ! Temperatura exacta
    real(dp)                 :: err_new                                             ! errores absolutos
    real(dp),    allocatable :: ux_old(:),ux_new(:),diag(:),diag_sup(:),diag_inf(:),err(:)
    20 format(E13.6,x,E13.6,x,E13.6);21 format(E13.6,x,E13.6)
    open(10,file='../results/result_02_implicit_vn.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat
    open(11,file='../results/result_02_implicit_vn_err.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(11file) = ',istat
    ! calculamos la parte entera del modulo y del resto del inverso de x_step_adim
    a=int(mod(1._dp,x_step_adim),sp);b=int(1._dp/x_step_adim);n=(a+b-1_sp)*(1_sp/(1_sp-a))
    ! cargamos diagonales central, superior e inferior
    allocate(diag(n+2),diag_sup(n+2),diag_inf(n+2))
    diag(:)=1._dp+2._dp*alpha;diag_inf(:)=-alpha;diag_sup(:)=-alpha
    ! cambiamos los extremos de las diagonales iferior y superior
    diag_inf(n+2)=diag_inf(n+2)*2._dp;diag_sup(1)=diag_sup(1)*2._dp
    diag_sup(n+2)=0._dp;diag_inf(1)=0._dp
    ! cargamos datos iniciales de temperatura
    allocate(ux_old(n+2))
    do i=1,n+2;ux_old(i)=cos(pi*real(i-1,dp)*x_step_adim);end do
    ! aplicamos método implicito, hacemos la evolución temporal (considerando los extremos)
    allocate(ux_new(n+2),err(t_write))
    err(:)=0._dp
    do j=1,(t_write-1)
        ux_old(1)=ux_old(1)+2._dp*alpha*x_step_adim*u0_Li
        ux_old(n+2)=ux_old(n+2)-2._dp*alpha*x_step_adim*u0_Lf
        call implicit_method(n+2,diag,diag_sup,diag_inf,ux_old(:),ux_new(:))
        ux_old(:)=ux_new(:)
        ! imprimimos errores absolutos para todo t
        do i=1,n+2
            call exact_solution(real(j,dp)*t_step_adim,real(i-1,dp)*x_step_adim,T_exact)
            err_new=abs(T_exact-ux_new(i));if (err_new>err(j)) err(j)=err_new
        end do;write(11,21) real(j,dp)*t_step_adim,err(j)
    end do
    ! escribimos valores luego de t_write pasos temporales
    call implicit_method(n+2,diag,diag_sup,diag_inf,ux_old(:),ux_new(:))
    do j=1,n+2
        call exact_solution(t_write_adim,real(j-1,dp)*x_step_adim,T_exact)
        err_new=abs(T_exact-ux_new(j));if (err_new>err(t_write)) err(t_write)=err_new
        write(10,20) real(j-1,dp)*x_step_adim,ux_new(j),abs((T_exact-ux_new(j))*(1._dp/T_exact))
    end do
    write(11,21) real(t_write,dp)*t_step_adim,err(t_write)
    close(10);deallocate(diag,diag_sup,diag_inf,ux_old,ux_new)
end program heateq_comparison_02
subroutine exact_solution(t_adim,x_adim,T_exact)
    use module_precision
    implicit none
    real(dp), intent(in)  :: t_adim,x_adim
    real(dp), intent(out) ::  T_exact
    real(dp), parameter   :: pi=4._dp*atan(1.0)
    T_exact=exp(-pi*pi*t_adim)*cos(pi*x_adim)
end subroutine exact_solution