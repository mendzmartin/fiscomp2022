! Comparación entre el método explícito y la solución exacta
program heateq_comparison_01
    use module_precision
    implicit none
    real(dp), allocatable   :: ux_old(:),ux_new(:)
    integer(dp), parameter  :: n=8_dp                   ! numero de elementos finitos
    real(dp), parameter     :: u0_Li=0._dp,u0_Lf=0._dp  ! condiciones de borde
    real(dp)                :: u0_ti=100._dp            ! condicion inicial
    real(dp)                :: param_t,param_x          ! parametros para adimensionalizar (espacial y temporal)
    real(dp)                :: D                        ! thermal diffusivity
    real(dp), parameter     :: alpha=1.E-7_dp           ! se asegura estabilidad (alpha < 0.5)
    real(dp)                :: T_exact                  ! temperatura exacta
    real(dp)                :: t_step_adim,x_step_adim  ! variables adimensionales (espacial y temporal)
    integer(dp)             :: i,j,k,istat              ! variables de loops
    real(dp), parameter     :: t_write_1_noadm=180._dp,t_write_2_noadm=1800._dp ! tiempos en dimensiones para escrituras
    integer(dp)             :: t_write_1,t_write_2      ! pasos temporales entre escrituras
    ! parametros para adimensionalizar
    D=237._dp*(1._dp/(900._dp*2700._dp))    ! D = K/C*rho [L^2]/[t]
    param_x=1._dp                         ! longitud caracteristica (long total barra en metros) [L]
    param_t=param_x*param_x*(1._dp/D)     ! tiempo caracteristico [t]
    x_step_adim=1._dp/(real(n,dp)+1._dp)
    t_step_adim=alpha*x_step_adim*x_step_adim
    write(*,*) 'D           = ',D
    write(*,*) 'param_t     = ',param_t
    write(*,*) 't_step_adim = ',t_step_adim
    write(*,*) 'x_step_adim = ',x_step_adim
    allocate(ux_old(n+2))
    ! inicializamos el vector espacial a tiempo inicial
    ux_old(1)=u0_Li
    ux_old(2:n+1)=u0_ti
    ux_old(n+2)=u0_Lf
    allocate(ux_new(n+2))
    ! ADIMENSIONALIZAMOS LOS TIEMPOS DE ESCRITURA
    t_write_1=int(t_write_1_noadm*(1._dp/(param_t*t_step_adim)),dp)
    t_write_2=int(t_write_2_noadm*(1._dp/(param_t*t_step_adim)),dp)
    write(*,*) t_write_1
    ! HACEMOS LA EVOLUCIÓN TEMPORAL (SIN CONSIDERAR LOS EXTREMOS)
    ux_new(1)=u0_Li
    ux_new(n+2)=u0_Lf
    do i=1,(t_write_1-1)
        do j=2,(n+1)
            ux_new(j)=(1._dp-2._dp*alpha)*ux_old(j)+alpha*(ux_old(j+1)+ux_old(j-1))
            ux_old(j)=ux_new(j)
        end do
    end do
    ! ESCRIBIMOS VALORES LUEGO DE t_write_1 PASOS TEMPORALES
    20 format (E9.2,x,E9.2)
    open(10,file='../results/result_02_aprox_explicit.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'Input/Output file. istat10 = ',istat
    open(11,file='../results/result_02_exact.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'Input/Output file. istat11 = ',istat
    !call exact_solution(0._dp,t_write_1_noadm*(1._dp/param_t),u0_Li,T_exact)
    call exact_solution(0._dp,t_write_1_noadm,D,u0_Li,T_exact)
    write(10,20) 0._dp,ux_new(1)
    write(11,20) 0._dp,T_exact
    do j=2,(n+1)
        ux_new(j)=(1._dp-2._dp*alpha)*ux_old(j)+alpha*(ux_old(j+1)+ux_old(j-1))
        ux_old(j)=ux_new(j)
        !call exact_solution(real(j-1,dp)*x_step_adim,t_write_1_noadm*(1._dp/param_t),u0_ti,T_exact)
        call exact_solution(real(j-1,dp)*x_step_adim,t_write_1_noadm,D,u0_ti,T_exact)
        write(10,20) real(j-1,dp)*x_step_adim,ux_new(j)
        write(11,20) real(j-1,dp)*x_step_adim,T_exact
    end do
    !call exact_solution(real(n+1,dp)*x_step_adim,t_write_1_noadm*(1._dp/param_t),u0_Lf,T_exact)
    call exact_solution(real(n+1,dp)*x_step_adim,t_write_1_noadm,D,u0_Lf,T_exact)
    write(10,20) real(n+1,dp)*x_step_adim,ux_new(n+2)
    write(11,20) real(n+1,dp)*x_step_adim,T_exact
    do i=t_write_1,(t_write_2-1)
        do j=2,(n+1)
            ux_new(j)=(1._dp-2._dp*alpha)*ux_old(j)+alpha*(ux_old(j+1)+ux_old(j-1))
            ux_old(j)=ux_new(j)
        end do
    end do
    ! ESCRIBIMOS VALORES LUEGO DE t_write_2 PASOS TEMPORALES
    !call exact_solution(0._dp,t_write_2_noadm*(1._dp/param_t),u0_Li,T_exact)
    call exact_solution(0._dp,t_write_2_noadm,D,u0_Li,T_exact)
    write(10,20) 0._dp,ux_new(1)
    write(11,20) 0._dp,T_exact
    do j=2,(n+1)
        ux_new(j)=(1._dp-2._dp*alpha)*ux_old(j)+alpha*(ux_old(j+1)+ux_old(j-1))
        ux_old(j)=ux_new(j)
        !call exact_solution(real(j-1,dp)*x_step_adim,t_write_2_noadm*(1._dp/param_t),u0_ti,T_exact)
        call exact_solution(real(j-1,dp)*x_step_adim,t_write_2_noadm,D,u0_ti,T_exact)
        write(10,20) real(j-1,dp)*x_step_adim,ux_new(j)
        write(11,20) real(j-1,dp)*x_step_adim,T_exact
    end do
    !call exact_solution(real(n+1,dp)*x_step_adim,t_write_2_noadm*(1._dp/param_t),u0_Lf,T_exact)
    call exact_solution(real(n+1,dp)*x_step_adim,t_write_2_noadm,D,u0_Lf,T_exact)
    write(10,20) real(n+1,dp)*x_step_adim,ux_new(n+2)
    write(11,20) real(n+1,dp)*x_step_adim,T_exact
    close(10)
    close(11)
    deallocate(ux_old,ux_new)
end program heateq_comparison_01
! subrutina para calcular solución exacta
! subroutine exact_solution(x_adm,t_adm,Temp0,Temp)
!     use module_precision
!     implicit none
!     ! variables de entrada/salida
!     real(dp), intent(in)    :: x_adm,t_adm,Temp0 ! ojo, x y t son adimensionales
!     real(dp), intent(out)   :: Temp
!     ! variables locales
!     real(dp),    parameter  :: pi=4._dp*atan(1._dp)
!     integer(dp), parameter  :: n_max=81_dp ! debe ser impar
!     real(dp)                :: kn
!     integer(dp)             :: n
!     Temp=0._dp
!     do n=1,n_max,2 ! recorro numeros impares
!         kn=real(n,dp)*pi
!         Temp=Temp+(1._dp/real(n,dp))*sin(kn*x_adm)*exp(-(kn*kn*t_adm))
!     end do
!     Temp=4._dp*Temp0*(1._dp/pi)*Temp
! end subroutine exact_solution
subroutine exact_solution(x,t,D,Temp0,Temp)
    use module_precision
    implicit none
    ! variables de entrada/salida
    real(dp), intent(in)    :: x,t,Temp0,D ! ojo, x y t no son adimensionales
    real(dp), intent(out)   :: Temp
    ! variables locales
    real(dp),    parameter  :: pi=4._dp*atan(1._dp)
    integer(dp), parameter  :: n_max=81_dp ! debe ser impar
    real(dp)                :: kn
    integer(dp)             :: n
    Temp=0._dp
    do n=1,n_max,2 ! recorro numeros impares
        kn=real(n,dp)*pi
        Temp=Temp+(1._dp/real(n,dp))*sin(kn*x)*exp(-(kn*kn*t*D))
    end do
    Temp=4._dp*Temp0*(1._dp/pi)*Temp
end subroutine exact_solution