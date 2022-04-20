program eq_calor
    use module_presition
    
    implicit none
    
    real(dp), allocatable   :: ux_old(:), ux_new(:)
    integer(sp), parameter  :: n=100_sp                    ! numero de elementos finitos
    real(dp), parameter     :: u0_Li=0._dp, u0_Lf=0._dp ! condiciones de borde
    real(dp)                :: u0_ti=100._dp            ! condicion inicial
    real(dp)                :: param_t, param_x
    real(dp)                :: D
    real(dp)                :: alpha
    real(dp)                :: t_step_adim, x_step_adim
    integer(sp)             :: i,j,k,istat
    integer(sp), parameter  :: t_write=300_sp                  ! pasos temporales entre escrituras
    real(dp)                :: factor

    ! parametros para adimensionalizar
    D=900._dp*2700._dp*(1._dp/237._dp)      ! thermal diffusivity, C*rho/K [t/L^2]
    param_x = 1._dp                         ! longitud caracteristica (long total barra) [L]
    write(*,'(A20,E10.4)') 'param_x = ', param_x
    param_t = param_x*param_x*D             ! tiempo caracteristico (paso temporal) [t]
    write(*,'(A20,E10.4)') 'param_t = ', param_t
    t_step_adim = 0.3_dp*(1._dp/param_t)
    x_step_adim=1._dp/(real(n,dp)+1._dp)

    alpha=t_step_adim*(1._dp/(x_step_adim*x_step_adim))
    ! verificamos condición de estabilidad debe ser menor a 1/2
    write(*,'(A20,E10.4)') 'alpha = ', alpha


    allocate(ux_old(n))
    ! inicializamos el vector espacial a tiempo inicial
    ux_old(1)=u0_Li
    ux_old(n)=u0_Lf
    do i=2,(n-1_sp)
        ux_old(i)=u0_ti
    end do

    allocate(ux_new(n))
    ! hacemos la evolución temporal (sin considerar los extremos)
    ! asignamos valores al vector espacial a tiempo t
    20 format (2(E10.2))
    open( 10, file = '../results/result_01.dat', status = 'replace', action = 'write', iostat = istat )
    if (istat == 0_sp ) then

        ! escribimos valores inicialmente
        write(10,20) 0._dp, ux_new(1)
        do j=2_sp,(n-1_sp)
            ux_new(j)=alpha*ux_old(j-1)+(1._dp-2._dp*alpha)*ux_old(j)+alpha*ux_old(j+1)
            write(10,20) real(j,dp)*x_step_adim, ux_new(j)
            ux_old(j)=ux_new(j)
        end do
        write(10,20) real(n,dp)*x_step_adim, ux_new(n)

        do k=1,50
            do i=1_sp,(t_write-1_sp)
                ux_new(1)=u0_Li
                do j=2_sp,(n-1_sp)
                    ux_new(j)=alpha*ux_old(j-1)+(1._dp-2._dp*alpha)*ux_old(j)+alpha*ux_old(j+1)
                    ux_old(j)=ux_new(j)
                end do
                ux_new(n)=u0_Lf
            end do

            ! escribimos valores luego de t_write pasos temporales
            write(10,20) 0._dp, ux_new(1)
            do j=2_sp,(n-1_sp)
                ux_new(j)=alpha*ux_old(j-1)+(1._dp-2._dp*alpha)*ux_old(j)+alpha*ux_old(j+1)
                write(10,20) real(j,dp)*x_step_adim, ux_new(j)
                ux_old(j)=ux_new(j)
            end do
            write(10,20) real(n,dp)*x_step_adim, ux_new(n)
        end do
    else
        write(*,*) 'Input/Output file. istat = ', istat
    end if
    close(10)

    ! controlar termalización en el centro de la barra
    write(*,'(A20,E10.4)') 'T(L/2,t_final) = ', ux_new(50)
    ! estimamos tamaño de archivo
    factor=10._dp ! este factor depende de (20 format (2(E10.2))), en este ejemplo factor es 10
    write(*,'(A30,E10.4,A10)') 'Tamaño aprox. de archivo = ', 2._dp*real(n,dp)*51._dp*factor, '[bytes]'

    deallocate(ux_old,ux_new)
end program eq_calor