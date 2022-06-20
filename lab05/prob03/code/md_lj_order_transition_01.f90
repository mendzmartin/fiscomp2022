! make clean && make md_lj_order_transition_01.o && ./md_lj_order_transition_01.o
program md_lj_order_transition_01
    use module_precision;use module_md_lennard_jones
    implicit none
    integer(sp), parameter   :: n_p=256_sp                             ! cantidad de partículasa
    real(dp),    parameter   :: delta_time=0.005_dp                    ! paso temporal
    integer(sp), parameter   :: time_eq=2000_sp,&                      ! pasos de equilibración
                                time_run=1000_sp                       ! pasos de evolucion en el estado estacionario
    real(dp),    parameter   :: T_adim_ref=1.0_dp                      ! temperatura de referencia adimensional
    real(dp),    parameter   :: r_cutoff=2.5_dp,mass=1._dp             ! radio de corte de interacciones y masa     
    real(dp),    allocatable :: x_vector(:),y_vector(:),z_vector(:)    ! componentes de las posiciones/particula
    real(dp),    allocatable :: vx_vector(:),vy_vector(:),vz_vector(:) ! componentes de la velocidad/particula
    real(dp),    allocatable :: force_x(:),force_y(:),force_z(:)       ! componentes de la fuerza/particula
    integer(sp)              :: i,j,k,index,istat                          ! loop index
    real(dp)                 :: T_adim                                 ! Temperatura
    real(dp)                 :: vx_mc,vy_mc,vz_mc                      ! componentes de la velocidad del centro de masas
    real(dp)                 :: time,time_end,time_start               ! tiempos de CPU
    logical                  :: pressure_switch,structure_factor_switch,&   ! variables para decidir escritura de datos
    diffusion_coeff_switch
    ! VARIABLES PARA REALIZAR BARRIDO DE DENSIDADES
    real(dp)                 :: density                                ! densidad (particulas/volumen)
    real(dp),    parameter   :: density_min=0.8_dp,density_max=1.2_dp  ! rango de densidades
    integer(sp), parameter   :: n_density=10_sp                        ! cantidad de densidades simuladas
    real(dp),    parameter   :: step_density=abs(density_max-density_min)*& ! paso de variació de densidades
                                             (1._dp/real(n_density-1,dp))
    ! VARIABLES PARA COMPUTAR PRESIÓN Y FACTOR DE ESTRUCTURA ESTÁTICO
    real(dp)                 :: press,press_med,var_press,err_press
    real(dp)                 :: s1_press,s2_press
    real(dp)                 :: Sk,Sk_med,var_Sk,err_Sk
    real(dp)                 :: s1_Sk,s2_Sk
    real(dp)                 :: D,D_med,var_D,err_D
    real(dp)                 :: s1_D,s2_D
    real(dp)                 :: msd_med,var_msd,err_msd
    real(dp)                 :: s1_msd,s2_msd
    ! VARIABLES PARA COMPUTAR COEFICIENTE DE DIFUSIÓN
    integer(sp), parameter   :: tau_max_corr=1000_sp                   ! pasos maximos de correlación
    real(dp),    allocatable :: wxx_matrix(:,:),wyy_matrix(:,:),&      ! matrices auxiliares para cálculo de msd
                                wzz_matrix(:,:)
    real(dp),    allocatable :: sum_wxx_vector(:),sum_wyy_vector(:),&  ! vectores auxiliares para cálculo de msd
                                sum_wzz_vector(:),counter_data(:)
    real(dp)                 :: msd                                    ! desplazamiento cuadrático medio
    integer(sp)              :: counter

    pressure_switch         =.false. ! escribir presión vs densidad
    structure_factor_switch =.false. ! escribir factor de estructura vs densidad
    diffusion_coeff_switch  =.true. ! escribir coeficiente de difusión vs densidad

    20 format(2(E12.4,x),x,E12.4);21 format(2(A12,x),x,A12)
    22 format(4(E12.4,x),x,E12.4);23 format(4(A12,x),x,A12)
    if (pressure_switch.eqv..true.) then
        open(10,file='../results/pressure_vs_density.dat',status='replace',action='write',iostat=istat)
        if (istat/=0) write(*,*) 'ERROR! istat(10file) = ',istat
        write(10,21) 'density','pressure','error'
    else if (structure_factor_switch.eqv..true.) then
        open(11,file='../results/struct_factor_vs_density.dat',status='replace',action='write',iostat=istat)
        if (istat/=0) write(*,*) 'ERROR! istat(11file) = ',istat
        write(11,21) 'density','S(k)_med','error'
    else if (diffusion_coeff_switch.eqv..true.) then
        open(12,file='../results/diffsuion_vs_density.dat',status='replace',action='write',iostat=istat)
        if (istat/=0) write(*,*) 'ERROR! istat(12file) = ',istat
        write(12,23) 'density','D','error','msd','error'
    end if

    call cpu_time(time_start)

    ! allocación de memoria
    allocate(x_vector(n_p),y_vector(n_p),z_vector(n_p))
    allocate(vx_vector(n_p),vy_vector(n_p),vz_vector(n_p))

    allocate(force_x(n_p),force_y(n_p),force_z(n_p))

    allocate(wxx_matrix(n_p,tau_max_corr),wyy_matrix(n_p,tau_max_corr),wzz_matrix(n_p,tau_max_corr))
    allocate(sum_wxx_vector(tau_max_corr),sum_wyy_vector(tau_max_corr),sum_wzz_vector(tau_max_corr))
    allocate(counter_data(tau_max_corr))

    do j=1,n_density
        ! seteo de variables y parámetros
        x_vector(:)=0._dp;y_vector(:)=0._dp;z_vector(:)=0._dp
        vx_vector(:)=0._dp;vy_vector(:)=0._dp;vz_vector(:)=0._dp
        force_x(:)=0._dp;force_y(:)=0._dp;force_z(:)=0._dp
        sum_wxx_vector(:)=0._dp;sum_wxx_vector(:)=0._dp;sum_wxx_vector(:)=0._dp
        counter=0;counter_data(:)=0_sp

        ! definimos r_cutoff en el rango [density_min;density_max]
        density=density_min+step_density*real(j-1,dp)

        ! generamos configuración inicial (FCC structure)
        call initial_lattice_configuration(n_p,density,x_vector,y_vector,z_vector,2)
        call mean_squared_displacement(n_p,x_vector,y_vector,z_vector,tau_max_corr,&
        wxx_matrix,wyy_matrix,wzz_matrix,sum_wxx_vector,sum_wyy_vector,sum_wzz_vector,&
        counter_data,counter)

        call md_initial_parameters(n_p,x_vector,y_vector,z_vector,&
        vx_vector,vy_vector,vz_vector,T_adim_ref,delta_time,density,mass)

        ! computamos fuerzas en el tiempo inicial
        call f_lj_total(x_vector,y_vector,z_vector,r_cutoff,n_p,density,force_x,force_y,force_z)

        ! TRANSITORIO
        index=0
        time=0._dp
        do i=1,time_eq
            index=index+1
            write(*,*) 'paso temporal =',index,' de',time_eq+time_run,j
            call rescaling_velocities(n_p,vx_vector,vy_vector,vz_vector,T_adim_ref,mass)
            call velocity_verlet(n_p,x_vector,y_vector,z_vector,&
            vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force_x,force_y,force_z)
            ! velocity center of mass to zero
            vx_mc=sum(vx_vector(:))*(1._dp/real(n_p,dp));vx_vector(:)=(vx_vector(:)-vx_mc)
            vy_mc=sum(vy_vector(:))*(1._dp/real(n_p,dp));vy_vector(:)=(vy_vector(:)-vy_mc)
            vz_mc=sum(vz_vector(:))*(1._dp/real(n_p,dp));vz_vector(:)=(vz_vector(:)-vz_mc)
            T_adim=temperature(n_p,mass,vx_vector,vy_vector,vz_vector)
            call mean_squared_displacement(n_p,x_vector,y_vector,z_vector,tau_max_corr,&
            wxx_matrix,wyy_matrix,wzz_matrix,sum_wxx_vector,sum_wyy_vector,sum_wzz_vector,&
            counter_data,counter)
            time=real(i,dp)*delta_time
        end do

        ! ESTACIONARIO
        press_med=0._dp;s1_press=0._dp;s2_press=0._dp
        Sk_med=0._dp;s1_Sk=0._dp;s2_Sk=0._dp

        time=0._dp
        do i=1,time_run
            index=index+1
            write(*,*) 'paso temporal =',index,' de',time_eq+time_run,j
            call rescaling_velocities(n_p,vx_vector,vy_vector,vz_vector,T_adim_ref,mass)
            call velocity_verlet(n_p,x_vector,y_vector,z_vector,&
            vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force_x,force_y,force_z)

            ! computamos observables,1er y 2do momento,valores medios y varianzas
            if (pressure_switch.eqv..true.) then
                press=pressure(n_p,density,mass,r_cutoff,x_vector,y_vector,z_vector,&
                vx_vector,vy_vector,vz_vector)
                s1_press=s1_press+press;s2_press=s2_press+press*press
                press_med=s1_press*(1._dp/real(i,dp))
                var_press=(real(i,dp)*s2_press-s1_press*s1_press)*(1._dp/real(i*i,dp))
            else if (structure_factor_switch.eqv..true.) then
                Sk=static_structure_factor(n_p,density,x_vector,y_vector,z_vector)
                s1_Sk=s1_Sk+Sk;s2_Sk=s2_Sk+Sk*Sk
                Sk_med=s1_Sk*(1._dp/real(i,dp))
                var_Sk=(real(i,dp)*s2_Sk-s1_Sk*s1_Sk)*(1._dp/real(i*i,dp))
            else if (diffusion_coeff_switch.eqv..true.) then
                call mean_squared_displacement(n_p,x_vector,y_vector,z_vector,tau_max_corr,&
                wxx_matrix,wyy_matrix,wzz_matrix,sum_wxx_vector,sum_wyy_vector,sum_wzz_vector,&
                counter_data,counter)
            end if     

            time=real(i,dp)*delta_time
        end do

        if (pressure_switch.eqv..true.) then
            ! computamos errores en el último paso
            err_press=(var_press*0.25_dp)*(1._dp/real(time_eq-1,dp))
            write(10,20) density,press_med,err_press
        else if (structure_factor_switch.eqv..true.) then
            ! computamos errores en el último paso
            err_Sk=(var_Sk*0.25_dp)*(1._dp/real(time_eq-1,dp))
            write(11,20) density,Sk_med,err_Sk
        else if (diffusion_coeff_switch.eqv..true.) then
            D_med=0._dp;s1_D=0._dp;s2_D=0._dp
            msd_med=0._dp;s1_msd=0._dp;s2_msd=0._dp
            do k=1,tau_max_corr
                time=real(k,dp)*delta_time
                ! computamos msd
                msd=(sum_wxx_vector(k)+sum_wyy_vector(k)+sum_wzz_vector(k))*&
                    (1._dp/real(counter_data(k),dp))*(1._dp/real(n_p,dp))
                ! computamos observables,1er y 2do momento,valores medios y varianzas
                D=msd*(1._dp/6._dp)*(1._dp/time)
                s1_D=s1_D+D;s2_D=s2_D+D*D
                D_med=s1_D*(1._dp/real(k,dp))
                var_D=(real(k,dp)*s2_D-s1_D*s1_D)*(1._dp/real(k*k,dp))

                s1_msd=s1_msd+msd;s2_msd=s2_msd+msd*msd
                msd_med=s1_msd*(1._dp/real(k,dp))
                var_msd=(real(k,dp)*s2_msd-s1_msd*s1_msd)*(1._dp/real(k*k,dp))
            end do
            ! computamos errores en el último paso
            err_D=(var_D*0.25_dp)*(1._dp/real(tau_max_corr-1,dp))
            err_msd=(var_msd*0.25_dp)*(1._dp/real(tau_max_corr-1,dp))
            write(12,22) density,D_med,err_D,msd_med,err_msd
        end if
    end do

    if (pressure_switch.eqv..true.) close(10)
    if (structure_factor_switch.eqv..true.) close(11)
    if (diffusion_coeff_switch.eqv..true.) close(12)

    deallocate(x_vector,y_vector,z_vector)
    deallocate(vx_vector,vy_vector,vz_vector)
    deallocate(force_x,force_y,force_z)
    deallocate(wxx_matrix,wyy_matrix,wzz_matrix)
    deallocate(sum_wxx_vector,sum_wyy_vector,sum_wzz_vector)
    deallocate(counter_data)

    call cpu_time(time_end)
    write(*,*) 'elapsed time = ',time_end-time_start,'[s]'
end program md_lj_order_transition_01

! subrutina para calcular el desplazamiento cuadrático medio
subroutine mean_squared_displacement(n_p,x_vector,y_vector,z_vector,tau_max_corr,&
    wxx_matrix,wyy_matrix,wzz_matrix,sum_wxx_vector,sum_wyy_vector,sum_wzz_vector,&
    counter_data,counter)
    use module_precision

    implicit none
    integer(sp), intent(in)    :: n_p                                       ! numero total de partículas
    integer(sp), intent(in)    :: tau_max_corr                              ! pasos maximos de autocorrelación
    real(dp),    intent(in)    :: x_vector(n_p),y_vector(n_p),z_vector(n_p) ! componentes del vector posición
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
        wxx_matrix(i,tau_corr_0)=x_vector(i)
        wyy_matrix(i,tau_corr_0)=y_vector(i)
        wzz_matrix(i,tau_corr_0)=z_vector(i)
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