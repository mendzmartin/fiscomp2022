program heateq_comparison_01
    use module_precision
    implicit none
    real(dp),    parameter   :: u0_Li=0._dp,u0_Lf=0._dp                             ! condiciones de contorno von neumann
    real(dp),    parameter   :: t_step_adim=0.001_dp,x_step_adim=0.05_dp            ! deltaT y deltaX adimensionales
    real(dp),    parameter   :: t_write_adim=1._dp                                  ! tiempo adimensional de evolución
    real(dp),    parameter   :: pi=4._dp*atan(1.0)
    real(dp),    parameter   :: D=237._dp*(1._dp/(900._dp*2700._dp))                ! thermal diffusivity(D=(K/C*rho)[L^2]/[t])
    real(dp),    parameter   :: param_x=1._dp                                       ! longitud caracteristica (long total barra en metros) [L]
    real(dp),    parameter   :: param_t=param_x*param_x*(1._dp/D)                   ! tiempo caracteristico (paso temporal en segundos) [t]
    real(dp),    parameter   :: alpha=t_step_adim*(1._dp/(x_step_adim*x_step_adim))
    integer(sp), parameter   :: t_write=int(t_write_adim*(1._dp/t_step_adim),sp)    ! pasos temporales de evolución
    integer(sp)              :: i,j,k,istat
    real(dp)                 :: T_exact                                             ! Temperatura exacta
    real(dp)                 :: err_new                                             ! errores absolutos
    integer(sp)              :: n,a,b                                               ! numero de elementos finitos y variables auxiliares
    real(dp),    allocatable :: ux_old(:),ux_new(:),err(:)
    20 format(E13.6,x,E13.6);21 format(E13.6,x,E13.6,x,E13.6)
    open(10,file='../results/result_02_explicit_vn.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat
    open(11,file='../results/result_02_exact_vn.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(11file) = ',istat
    open(12,file='../results/result_02_explicit_vn_err.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(12file) = ',istat
    ! verificamos condición de estabilidad debe ser menor a 1/2
    write(*,'(A20,E10.4)') 'alpha = ',alpha; write(*,'(A20,E10.4)') 'param_t = ',param_t
    ! calculamos la parte entera del modulo y del resto del inverso de x_step_adim
    a=int(mod(1._dp,x_step_adim),sp);b=int(1._dp/x_step_adim);n=(a+b-1_sp)*(1_sp/(1_sp-a))
    allocate(ux_old(n+2))
    ! inicializamos el vector espacial a tiempo inicial
    do i=1,n+2;ux_old(i)=cos(pi*real(i-1,dp)*x_step_adim);end do
    ! hacemos la evolución temporal (considerando los extremos)
    allocate(ux_new(n+2),err(t_write))
    err(:)=0._dp
    do i=1,(t_write-1)
        ux_new(1)=(1._dp-2._dp*alpha)*ux_old(1)+2._dp*alpha*(ux_old(2)-x_step_adim*u0_Li)
        ! imprimimos errores absolutos para todo t
        call exact_solution(real(i,dp)*t_step_adim,0._dp,T_exact)
        err_new=abs(T_exact-ux_new(1));if (err_new>err(i)) err(i)=err_new
        do j=2,(n+1)
            ux_new(j)=(1._dp-2._dp*alpha)*ux_old(j)+alpha*(ux_old(j+1)+ux_old(j-1))
            call exact_solution(real(i,dp)*t_step_adim,real(j-1,dp)*x_step_adim,T_exact)
            err_new=abs(T_exact-ux_new(j));if (err_new>err(i)) err(i)=err_new
        end do
        ux_new(n+2)=(1._dp-2._dp*alpha)*ux_old(n+2)+2._dp*alpha*(ux_old(n+1)+x_step_adim*u0_Lf)
        ux_old(1:n+2)=ux_new(1:n+2)
        call exact_solution(real(i,dp)*t_step_adim,real(n+1,dp)*x_step_adim,T_exact)
        err_new=abs(T_exact-ux_new(j));if (err_new>err(i)) err(i)=err_new
        write(12,21) real(i,dp)*t_step_adim,err(i)
    end do
    ! escribimos valores luego de t_write pasos temporales
    ux_new(1)=(1._dp-2._dp*alpha)*ux_old(1)+2._dp*alpha*(ux_old(2)-x_step_adim*u0_Li)
    call exact_solution(t_write_adim,0._dp,T_exact)
    err_new=abs(T_exact-ux_new(1));if (err_new>err(t_write)) err(t_write)=err_new
    write(10,21) 0._dp,ux_new(1),abs((T_exact-ux_new(1))*(1._dp/T_exact));write(11,20) 0._dp,T_exact
    do j=2,(n+1)
        ux_new(j)=(1._dp-2._dp*alpha)*ux_old(j)+alpha*(ux_old(j+1)+ux_old(j-1))
        call exact_solution(t_write_adim,real(j-1,dp)*x_step_adim,T_exact)
        err_new=abs(T_exact-ux_new(j));if (err_new>err(t_write)) err(t_write)=err_new
        write(10,21) real(j-1,dp)*x_step_adim,ux_new(j),abs((T_exact-ux_new(j))*(1._dp/T_exact))
        write(11,20) real(j-1,dp)*x_step_adim,T_exact
    end do
    ux_new(n+2)=(1._dp-2._dp*alpha)*ux_old(n+2)+2._dp*alpha*(ux_old(n+1)+x_step_adim*u0_Lf)
    call exact_solution(t_write_adim,real(n+1,dp)*x_step_adim,T_exact)
    err_new=abs(T_exact-ux_new(j));if (err_new>err(t_write)) err(t_write)=err_new
    write(10,21) real(n+1,dp)*x_step_adim,ux_new(n+2),abs((T_exact-ux_new(n+2))*(1._dp/T_exact))
    write(11,20) real(n+1,dp)*x_step_adim,T_exact;write(12,20) t_write_adim,err(t_write)
    close(10);deallocate(ux_old,ux_new)
end program heateq_comparison_01
subroutine exact_solution(t_adim,x_adim,T_exact)
    use module_precision
    implicit none
    real(dp), intent(in)  :: t_adim,x_adim
    real(dp), intent(out) ::  T_exact
    real(dp), parameter   :: pi=4._dp*atan(1.0)
    T_exact=exp(-pi*pi*t_adim)*cos(pi*x_adim)
end subroutine exact_solution