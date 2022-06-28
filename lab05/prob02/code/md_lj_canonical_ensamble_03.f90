! make clean && make md_lj_canonical_ensamble_03.o && ./md_lj_canonical_ensamble_03.o
program md_lj_canonical_ensamble_03
    use module_precision;use module_md_lennard_jones
    implicit none
    integer(sp), parameter   :: n_p=500_sp                             ! cantidad de partículasa
    real(dp),    parameter   :: delta_time=0.005_dp                    ! paso temporal
    integer(sp), parameter   :: time_eq=2000_sp,&                      ! pasos de equilibración
                                time_run=15000_sp                       ! pasos de evolucion en el estado estacionario
    integer(sp), parameter   :: tau_max_corr=1000_sp                   ! pasos maximos de correlación
    real(dp),    parameter   :: T_adim_ref=1.0_dp                      ! temperatura de referencia adimensional
    real(dp),    parameter   :: r_cutoff=2.5_dp,mass=1._dp             ! radio de corte de interacciones y masa     
    real(dp),    allocatable :: x_vector(:),y_vector(:),z_vector(:)    ! componentes de las posiciones/particula
    real(dp),    allocatable :: x_vector_noPBC(:),y_vector_noPBC(:),&  ! componentes de la posición sin PBC
                                z_vector_noPBC(:)
    real(dp),    allocatable :: vx_vector(:),vy_vector(:),vz_vector(:) ! componentes de la velocidad/particula
    real(dp),    allocatable :: force_x(:),force_y(:),force_z(:)       ! componentes de la fuerza/particula
    real(dp),    allocatable :: wxx_matrix(:,:),wyy_matrix(:,:),&      ! matrices auxiliares para cálculo de msd
                                wzz_matrix(:,:)
    real(dp),    allocatable :: sum_wxx_vector(:),sum_wyy_vector(:),&  ! vectores auxiliares para cálculo de msd
                                sum_wzz_vector(:),counter_data(:)
    integer(sp)              :: i,j,index,istat,counter                 ! loop index
    real(dp)                 :: T_adim                                 ! Temperatura
    real(dp)                 :: vx_mc,vy_mc,vz_mc                      ! componentes de la velocidad del centro de masas
    real(dp)                 :: time,time_end,time_start               ! tiempos de CPU
    real(dp)                 :: msd                                    ! desplazamiento cuadrático medio
    !real(dp),    parameter   :: density=0.8_dp                         ! densidad (particulas/volumen)
    real(dp),    parameter   :: density=1.2_dp                         ! densidad (particulas/volumen)

    ! DESCOMENTAR PARA density=0.8
    !open(10,file='../results/msd_rho1.dat',status='replace',action='write',iostat=istat)
    ! DESCOMENTAR PARA density=1.2
    open(10,file='../results/msd_rho2.dat',status='replace',action='write',iostat=istat)
    if (istat/=0) write(*,*) 'ERROR! istat(11file) = ',istat
    24 format(E12.4,x,E12.4);25 format(A12,x,A12)
    write(10,25) 'time','msd'

    call cpu_time(time_start)

    allocate(x_vector(n_p),y_vector(n_p),z_vector(n_p))
    x_vector(:)=0._dp;y_vector(:)=0._dp;z_vector(:)=0._dp
    allocate(x_vector_noPBC(n_p),y_vector_noPBC(n_p),z_vector_noPBC(n_p))
    x_vector_noPBC(:)=0._dp;y_vector_noPBC(:)=0._dp;z_vector_noPBC(:)=0._dp

    ! generamos configuración inicial (FCC structure)
    call initial_lattice_configuration(n_p,density,x_vector,y_vector,z_vector,2)
    x_vector_noPBC(:)=x_vector(:);y_vector_noPBC(:)=y_vector(:);z_vector_noPBC(:)=z_vector(:)

    counter=0_sp
    allocate(wxx_matrix(n_p,tau_max_corr),wyy_matrix(n_p,tau_max_corr),wzz_matrix(n_p,tau_max_corr))
    allocate(sum_wxx_vector(tau_max_corr),sum_wyy_vector(tau_max_corr),sum_wzz_vector(tau_max_corr))
    allocate(counter_data(tau_max_corr))

    sum_wxx_vector(:)=0._dp;sum_wxx_vector(:)=0._dp;sum_wxx_vector(:)=0._dp
    counter_data(:)=0._dp

    call mean_squared_displacement(n_p,x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,tau_max_corr,&
    wxx_matrix,wyy_matrix,wzz_matrix,sum_wxx_vector,sum_wyy_vector,sum_wzz_vector,&
    counter_data,counter)

    allocate(vx_vector(n_p),vy_vector(n_p),vz_vector(n_p))
    vx_vector(:)=0._dp;vy_vector(:)=0._dp;vz_vector(:)=0._dp
    call md_initial_parameters(n_p,x_vector,y_vector,z_vector,&
    vx_vector,vy_vector,vz_vector,T_adim_ref,delta_time,density,mass)

    ! computamos fuerzas en el tiempo inicial
    allocate(force_x(n_p),force_y(n_p),force_z(n_p))
    force_x(:)=0._dp;force_y(:)=0._dp;force_z(:)=0._dp
    call f_lj_total(x_vector,y_vector,z_vector,r_cutoff,n_p,density,force_x,force_y,force_z)

    ! TRANSITORIO
    index=0
    time=0._dp
    do i=1,time_eq
        index=index+1
        write(*,*) 'paso temporal =',index,' de',time_eq+time_run
        call rescaling_velocities(n_p,vx_vector,vy_vector,vz_vector,T_adim_ref,mass)
        call velocity_verlet(n_p,x_vector,y_vector,z_vector,&
        x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
        vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force_x,force_y,force_z)
        ! velocity center of mass to zero
        vx_mc=sum(vx_vector(:))*(1._dp/real(n_p,dp));vx_vector(:)=(vx_vector(:)-vx_mc)
        vy_mc=sum(vy_vector(:))*(1._dp/real(n_p,dp));vy_vector(:)=(vy_vector(:)-vy_mc)
        vz_mc=sum(vz_vector(:))*(1._dp/real(n_p,dp));vz_vector(:)=(vz_vector(:)-vz_mc)
        T_adim=temperature(n_p,mass,vx_vector,vy_vector,vz_vector)
        time=real(i,dp)*delta_time
        ! call mean_squared_displacement(n_p,x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,tau_max_corr,&
        ! wxx_matrix,wyy_matrix,wzz_matrix,sum_wxx_vector,sum_wyy_vector,sum_wzz_vector,&
        ! counter_data,counter)
    end do

    ! ESTACIONARIO
    time=0._dp
    do i=1,time_run
        index=index+1
        write(*,*) 'paso temporal =',index,' de',time_eq+time_run
        call rescaling_velocities(n_p,vx_vector,vy_vector,vz_vector,T_adim_ref,mass)
        call velocity_verlet(n_p,x_vector,y_vector,z_vector,&
        x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
        vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force_x,force_y,force_z)
        time=(real(time_eq,dp)+real(i,dp))*delta_time
        call mean_squared_displacement(n_p,x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,tau_max_corr,&
        wxx_matrix,wyy_matrix,wzz_matrix,sum_wxx_vector,sum_wyy_vector,sum_wzz_vector,&
        counter_data,counter)
    end do

    do j=1,tau_max_corr
        time=real(j,dp)*delta_time
        msd=(sum_wxx_vector(j)+sum_wyy_vector(j)+sum_wzz_vector(j))*(1._dp/real(counter_data(j),dp))*(1._dp/real(n_p,dp))
        write(10,24) time,msd
    end do
    close(10)

    deallocate(x_vector,y_vector,z_vector)
    deallocate(vx_vector,vy_vector,vz_vector)
    deallocate(force_x,force_y,force_z)

    deallocate(wxx_matrix,wyy_matrix,wzz_matrix)
    deallocate(sum_wxx_vector,sum_wyy_vector,sum_wzz_vector)
    deallocate(counter_data)

    call cpu_time(time_end)
    write(*,*) 'elapsed time = ',time_end-time_start,'[s]'
end program md_lj_canonical_ensamble_03

! subrutina para calcular el desplazamiento cuadrático medio
subroutine mean_squared_displacement(n_p,x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,tau_max_corr,&
    wxx_matrix,wyy_matrix,wzz_matrix,sum_wxx_vector,sum_wyy_vector,sum_wzz_vector,&
    counter_data,counter)
    use module_precision

    implicit none
    integer(sp), intent(in)    :: n_p                                       ! numero total de partículas
    integer(sp), intent(in)    :: tau_max_corr                              ! pasos maximos de autocorrelación
    real(dp),    intent(in)    :: x_vector_noPBC(n_p),y_vector_noPBC(n_p),&
                                  z_vector_noPBC(n_p) ! componentes del vector posición
    real(dp),    intent(inout) :: wxx_matrix(n_p,tau_max_corr),&            ! matrices de acumulación
                                  wyy_matrix(n_p,tau_max_corr),&
                                  wzz_matrix(n_p,tau_max_corr)
    real(dp),    intent(inout) :: sum_wxx_vector(tau_max_corr),&            ! vectores de sumas auxiliares
                                  sum_wyy_vector(tau_max_corr),&
                                  sum_wzz_vector(tau_max_corr)
    real(dp),    intent(inout) :: counter_data(tau_max_corr)                ! contador de datos
    integer(sp), intent(inout) :: counter                                   ! contador de entradas
    
    integer(sp)            :: i,j
    integer(sp)            :: tau_corr_0,tau_corr_t ! tiempos de correlación
    integer(sp), parameter :: nmax_tau_corr_0=10_sp ! maximo número de tau_corr_0 que almacenamos

    counter=counter+1                        ! numero de veces que entro a la subrutina
    tau_corr_0=mod(counter-1,tau_max_corr)+1 ! tiempo de correlación actual tau_corr_0={1,2,...,tau_max_corr}

    ! guardamos cíclicamente los últimos tau_max_corr valores
    !  de las componentes x,y,z de cada partícula
    do i=1,n_p
        wxx_matrix(i,tau_corr_0)=x_vector_noPBC(i)
        wyy_matrix(i,tau_corr_0)=y_vector_noPBC(i)
        wzz_matrix(i,tau_corr_0)=z_vector_noPBC(i)
    end do

    if ((mod(counter,nmax_tau_corr_0)==0).and.(counter>tau_max_corr)) then
        do j=1,tau_max_corr
            tau_corr_t=mod(counter-j,tau_max_corr)+1
            do i=1,n_p
                sum_wxx_vector(j)=sum_wxx_vector(j)+(wxx_matrix(i,tau_corr_0)-wxx_matrix(i,tau_corr_t))*&
                                                    (wxx_matrix(i,tau_corr_0)-wxx_matrix(i,tau_corr_t))
                sum_wyy_vector(j)=sum_wyy_vector(j)+(wyy_matrix(i,tau_corr_0)-wyy_matrix(i,tau_corr_t))*&
                                                    (wyy_matrix(i,tau_corr_0)-wyy_matrix(i,tau_corr_t))
                sum_wzz_vector(j)=sum_wzz_vector(j)+(wzz_matrix(i,tau_corr_0)-wzz_matrix(i,tau_corr_t))*&
                                                    (wzz_matrix(i,tau_corr_0)-wzz_matrix(i,tau_corr_t))
            end do
            ! actualizamos el contador de datos para cada tiempo de correlación
            counter_data(j)=counter_data(j)+1._dp
        end do
    end if
end subroutine mean_squared_displacement