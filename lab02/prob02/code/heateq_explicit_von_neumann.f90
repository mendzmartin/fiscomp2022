program heateq_explicit_von_neumann
    use module_precision
    implicit none
    real(dp),    allocatable :: ux_old(:),ux_new(:)
    integer(sp), parameter   :: n=19_sp                  ! numero de elementos finitos
    real(dp),    parameter   :: u0_Li=0._dp,u0_Lf=0._dp  ! condiciones de contorno von neumann
    real(dp)                 :: pi=4._dp*atan(1.0_dp)
    real(dp)                 :: D
    real(dp)                 :: param_t,param_x          ! parametros para adimensionalizar (espacial y temporal)
    real(dp)                 :: alpha
    real(dp)                 :: t_step_adim,x_step_adim
    integer(sp)              :: i,j,k,istat
    integer(sp), parameter   :: t_write=2_sp                  ! pasos temporales entre escrituras
    20 format(E13.6,x,E13.6)
    open(10,file='../results/result_01_explicit_vn.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat
    ! ANALISIS DIMENSIONAL
    D=237._dp*(1._dp/(900._dp*2700._dp))    ! D = K/C*rho [L^2]/[t]
    param_x=1._dp                         ! longitud caracteristica (long total barra en metros) [L]
    param_t=param_x*param_x*(1._dp/D)     ! tiempo caracteristico [t]
    t_step_adim=0.001
    x_step_adim=1._dp/(real(n,dp)+1._dp)
    alpha=t_step_adim*(1._dp/(x_step_adim*x_step_adim))
    ! verificamos condición de estabilidad debe ser menor a 1/2
    write(*,'(A20,E10.4)') 'alpha = ',alpha
    write(*,*) 'param_t = ', param_t
    allocate(ux_old(n+2))
    ! inicializamos el vector espacial a tiempo inicial
    do i=1,n+2
        ux_old(i)=cos(pi*real(i-1,dp)*x_step_adim)
        write(10,20) real(i-1,dp)*x_step_adim,ux_old(i)
    end do
    allocate(ux_new(n+2))
    ! hacemos la evolución temporal (sin considerar los extremos)
    ! asignamos valores al vector espacial a tiempo t
    do k=1,10 ! evolucionamos (50*t_write) pasos temporales
        do i=1,(t_write-1)
            ux_new(1)=(1._dp-2._dp*alpha)*ux_old(1)+2._dp*alpha*(ux_old(2)-x_step_adim*u0_Li)
            do j=2,(n+1)
                ux_new(j)=(1._dp-2._dp*alpha)*ux_old(j)+alpha*(ux_old(j+1)+ux_old(j-1))
            end do
            ux_new(n+2)=(1._dp-2._dp*alpha)*ux_old(n+2)+2._dp*alpha*(ux_old(n+1)+x_step_adim*u0_Lf)
            ux_old(1:n+2)=ux_new(1:n+2)
        end do
        ! escribimos valores luego de t_write pasos temporales
        ux_new(1)=(1._dp-2._dp*alpha)*ux_old(1)+2._dp*alpha*(ux_old(2)-x_step_adim*u0_Li)
        write(10,20) 0._dp,ux_new(1)
        do j=2,(n+1)
            ux_new(j)=(1._dp-2._dp*alpha)*ux_old(j)+alpha*(ux_old(j+1)+ux_old(j-1))
            write(10,20) real(j-1,dp)*x_step_adim,ux_new(j)
        end do
        ux_new(n+2)=(1._dp-2._dp*alpha)*ux_old(n+2)+2._dp*alpha*(ux_old(n+1)+x_step_adim*u0_Lf)
        write(10,20) real(n+1,dp)*x_step_adim,ux_new(n+2)
    end do
    close(10)
    ! controlar termalización en el centro de la barra
    !write(*,'(A20,E10.4)') 'T(L/2,t_final) = ', ux_new(50)
    deallocate(ux_old,ux_new)
end program heateq_explicit_von_neumann