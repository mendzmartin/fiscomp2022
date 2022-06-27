! make clean && make brownian_dynamic_lennard_jones_02.o && ./brownian_dynamic_lennard_jones_02.o
program brownian_dynamic_lennard_jones_02
    use module_precision;use module_bd_lennard_jones
    implicit none
    ! VARIABLES y PARAMETROS GENERALES
    integer(sp), parameter   :: n_p=256_sp                             ! cantidad de partículasa
    real(dp),    parameter   :: delta_time=0.001_dp                    ! paso temporal
    integer(sp), parameter   :: time_eq=15000_sp,&                        ! pasos de equilibración
                                time_run=15000_sp,&                    ! pasos de evolucion en el estado estacionario
                                ensamble_step=50_sp                    ! pasos de evolución para promedio en ensamble
    real(dp),    parameter   :: T_adim_ref=1.1_dp                       ! temperatura de referencia adimensional
    real(dp),    parameter   :: r_cutoff=2.5_dp,mass=1._dp             ! radio de corte de interacciones y masa    
    real(dp),    parameter   :: dinamic_viscosity=2.87_dp
    real(dp),    parameter   :: pi=4._dp*atan(1._dp)
    real(dp),    parameter   :: diffusion_coeff=T_adim_ref*&
                                 (1._dp/(3._dp*pi*dinamic_viscosity)) 
    real(dp),    allocatable :: x_vector(:),y_vector(:),z_vector(:)    ! componentes de las posiciones/particula
    real(dp),    allocatable :: force_x(:),force_y(:),force_z(:)       ! componentes de la fuerza/particula
    integer(sp)              :: i,j,k,istat,index                        ! loop index
    real(dp)                 :: time_end,time_start                    ! tiempos de CPU
    ! VARIABLES PARA REALIZAR BARRIDO DE DENSIDADES
    real(dp)                 :: density                                ! densidad (particulas/volumen)
    real(dp),    parameter   :: density_min=0.2_dp,density_max=1.05_dp ! rango de densidades
    integer(sp), parameter   :: n_density=20_sp                        ! cantidad de densidades simuladas
    real(dp),    parameter   :: step_density=abs(density_max-density_min)*& ! paso de variació de densidades
                                             (1._dp/real(n_density-1,dp))
    ! VARIABLES LOGICAS PARA DECIDIR QUÉ ESCRIBIR
    logical, parameter       :: energies_switch=.true.,&    ! escribir energías en el estado estacionario
                                msd_switch=.true.         ! escribir coeficiente de difusión vs densidad
    ! VARIABLES PARA COMPUTAR ENERGÍA POTENCIAL
    real(dp)                 :: U_adim,time,press  ! observables
    real(dp)                 :: U_med,var_U,err_U
    real(dp)                 :: s1_U,s2_U
    ! VARIABLES PARA COMPUTAR PRESIÓN OSMÓTICA
    real(dp)                 :: press_med,var_press,err_press
    real(dp)                 :: s1_press,s2_press
    ! VARIABLES PARA COMPUTAR DESPLAZAMIENTO CUADRÁTICO MEDIO
    integer(sp), parameter   :: tau_max_corr=100_sp                   ! pasos maximos de correlación
    real(dp),    allocatable :: x_vector_noPBC(:),y_vector_noPBC(:),&  ! componentes de la posición sin PBC
                                z_vector_noPBC(:)
    real(dp),    allocatable :: wxx_matrix(:,:),wyy_matrix(:,:),&      ! matrices auxiliares para cálculo de msd
                                wzz_matrix(:,:)
    real(dp),    allocatable :: sum_wxx_vector(:),sum_wyy_vector(:),&  ! vectores auxiliares para cálculo de msd
                                sum_wzz_vector(:)
    integer(sp), allocatable :: counter_data(:)
    real(dp)                 :: msd,msd_med,var_msd,err_msd            ! desplazamiento cuadrático medio
    real(dp)                 :: s1_msd,s2_msd
    integer(sp)              :: counter
    ! VARIABLES PARA COMPUTAR COEFICIENTE DE DIFUSIÓN
    real(dp)                 :: D,D_med,var_D,err_D
    real(dp)                 :: s1_D,s2_D

    call cpu_time(time_start)

    allocate(x_vector(n_p),y_vector(n_p),z_vector(n_p))
    allocate(x_vector_noPBC(n_p),y_vector_noPBC(n_p),z_vector_noPBC(n_p))
    allocate(force_x(n_p),force_y(n_p),force_z(n_p))
    allocate(wxx_matrix(n_p,tau_max_corr),wyy_matrix(n_p,tau_max_corr),wzz_matrix(n_p,tau_max_corr))
    allocate(sum_wxx_vector(tau_max_corr),sum_wyy_vector(tau_max_corr),sum_wzz_vector(tau_max_corr))
    allocate(counter_data(tau_max_corr))

    ! ESCRIBIMOS DATOS
    if (energies_switch.eqv..true.) then
        open(12,file='../results/energies_vs_density.dat',status='replace',action='write',iostat=istat)
        write(*,*) 'istat(12file) = ',istat;write(12,"(4(A12,x),A12)") 'density','U_med','err_U','press_med','err_press'
    end if
    if (msd_switch.eqv..true.) then
        open(13,file='../results/msd_vs_density.dat',status='replace',action='write',iostat=istat)
        if (istat/=0) write(*,*) 'ERROR! istat(12file) = ',istat
        write(13,"(4(A12,x),A12)") 'density','msd_med','err_msd','D/D0_med','err_D/D0'
    end if

    do k=1,n_density
        print*, k

        ! seteo de variables y parámetros
        x_vector(:)=0._dp;y_vector(:)=0._dp;z_vector(:)=0._dp
        x_vector_noPBC(:)=0._dp;y_vector_noPBC(:)=0._dp;z_vector_noPBC(:)=0._dp
        force_x(:)=0._dp;force_y(:)=0._dp;force_z(:)=0._dp
        sum_wxx_vector(:)=0._dp;sum_wxx_vector(:)=0._dp;sum_wxx_vector(:)=0._dp
        counter=0_sp;counter_data(:)=0_sp

        ! definimos denisty en el rango [density_min;density_max]
        density=density_min+step_density*real(k-1,dp)

        ! generamos configuración inicial (FCC structure)
        call initial_lattice_configuration(n_p,density,x_vector,y_vector,z_vector,2)
        ! computamos fuerzas en el tiempo inicial
        call f_lj_total(x_vector,y_vector,z_vector,r_cutoff,n_p,density,force_x,force_y,force_z)

        if (msd_switch.eqv..true.) then
            call mean_squared_displacement(n_p,x_vector,y_vector,z_vector,tau_max_corr,&
                wxx_matrix,wyy_matrix,wzz_matrix,sum_wxx_vector,sum_wyy_vector,sum_wzz_vector,&
                counter_data,counter)
        end if

        ! TRANSITORIO
        do i=1,time_eq
            call evolution_bd(n_p,x_vector,y_vector,z_vector,&
                x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
                delta_time,mass,r_cutoff,density,force_x,force_y,force_z,&
                dinamic_viscosity,diffusion_coeff)
            if ((mod(i,ensamble_step)==0_sp).and.(msd_switch.eqv..true.)) then
                call mean_squared_displacement(n_p,x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,tau_max_corr,&
                    wxx_matrix,wyy_matrix,wzz_matrix,sum_wxx_vector,sum_wyy_vector,sum_wzz_vector,&
                    counter_data,counter)
            end if
        end do

        ! ESTACIONARIO
        if (energies_switch.eqv..true.) then
            U_med=0._dp;s1_U=0._dp;s2_U=0._dp
            press_med=0._dp;s1_press=0._dp;s2_press=0._dp
        end if

        ! seteamos variables
        time=0._dp;i=0_sp

        do j=1,time_run
            call evolution_bd(n_p,x_vector,y_vector,z_vector,&
                x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
                delta_time,mass,r_cutoff,density,force_x,force_y,force_z,&
                dinamic_viscosity,diffusion_coeff)

            ! ESCRIBIMOS DATOS
            if (mod(j,ensamble_step)==0_sp) then
                i=i+1_sp
                time=real(i,dp)*delta_time
                if (energies_switch.eqv..true.) then
                    U_adim=u_lj_total(n_p,x_vector,y_vector,z_vector,r_cutoff,density)
                    U_adim=U_adim*(1._dp/real(n_p,dp))
                    press=osmotic_pressure(n_p,density,mass,r_cutoff,x_vector,y_vector,z_vector)

                    s1_U=s1_U+U_adim;s2_U=s2_U+U_adim*U_adim
                    U_med=s1_U*(1._dp/real(i,dp))
                    var_U=(real(i,dp)*s2_U-s1_U*s1_U)*(1._dp/real(i*i,dp))

                    s1_press=s1_press+press;s2_press=s2_press+press*press
                    press_med=s1_press*(1._dp/real(i,dp))
                    var_press=(real(i,dp)*s2_press-s1_press*s1_press)*(1._dp/real(i*i,dp))
                end if
                if (msd_switch.eqv..true.) then
                    call mean_squared_displacement(n_p,x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,tau_max_corr,&
                    wxx_matrix,wyy_matrix,wzz_matrix,sum_wxx_vector,sum_wyy_vector,sum_wzz_vector,&
                    counter_data,counter)
                end if
            end if
        end do

        if (energies_switch.eqv..true.) then
            ! computamos errores en el último paso
            !err_U=(var_U*0.25_dp)*(1._dp/real(i-1,dp))
            !err_press=(var_press*0.25_dp)*(1._dp/real(i-1,dp))
            err_U=sqrt(var_U*(1._dp/real(i-1,dp)))
            err_press=sqrt(var_press*(1._dp/real(i-1,dp)))

            write(12,"(4(E12.4,x),E12.4)") density,U_med,err_U,press_med,err_press
        end if

        if (msd_switch.eqv..true.) then
            msd_med=0._dp;s1_msd=0._dp;s2_msd=0._dp
            D_med=0._dp;s1_D=0._dp;s2_D=0._dp
            do i=1,tau_max_corr
                if (counter_data(i)==0_sp) stop
                time=real(i,dp)*delta_time
                ! computamos msd
                msd=(sum_wxx_vector(i)+sum_wyy_vector(i)+sum_wzz_vector(i))*&
                    (1._dp/real(counter_data(i),dp))*(1._dp/real(n_p,dp))
                D=msd*(1._dp/(6._dp*time))*(1._dp/diffusion_coeff)
                ! computamos observables,1er y 2do momento,valores medios y varianzas
                s1_msd=s1_msd+msd;s2_msd=s2_msd+msd*msd
                msd_med=s1_msd*(1._dp/real(i,dp))
                var_msd=(real(i,dp)*s2_msd-s1_msd*s1_msd)*(1._dp/real(i*i,dp))

                s1_D=s1_D+D;s2_D=s2_D+D*D
                D_med=s1_D*(1._dp/real(i,dp))
                var_D=(real(i,dp)*s2_D-s1_D*s1_D)*(1._dp/real(i*i,dp))
            end do
            ! computamos errores en el último paso
            !err_msd=(var_msd*0.25_dp)*(1._dp/real(tau_max_corr-1,dp))
            err_msd=sqrt(var_msd*(1._dp/real(tau_max_corr-1,dp)))
            err_D=sqrt(var_D*(1._dp/real(tau_max_corr-1,dp)))
            write(13,"(4(E12.4,x),E12.4)") density,msd_med,err_msd,D_med,err_D
        end if
    end do

    if (energies_switch.eqv..true.) close(12)
    if (msd_switch.eqv..true.) close(13)

    deallocate(x_vector,y_vector,z_vector)
    deallocate(x_vector_noPBC,y_vector_noPBC,z_vector_noPBC)
    deallocate(force_x,force_y,force_z)
    deallocate(wxx_matrix,wyy_matrix,wzz_matrix)
    deallocate(sum_wxx_vector,sum_wyy_vector,sum_wzz_vector)
    deallocate(counter_data)

    call cpu_time(time_end)
    write(*,*) 'elapsed time = ',time_end-time_start,'[s]'

end program brownian_dynamic_lennard_jones_02

! subrutina para calcular el desplazamiento cuadrático medio
subroutine mean_squared_displacement(n_p,x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,tau_max_corr,&
    wxx_matrix,wyy_matrix,wzz_matrix,sum_wxx_vector,sum_wyy_vector,sum_wzz_vector,&
    counter_data,counter)
    use module_precision

    implicit none
    integer(sp), intent(in)    :: n_p                                       ! numero total de partículas
    integer(sp), intent(in)    :: tau_max_corr                              ! pasos maximos de autocorrelación
    real(dp),    intent(in)    :: x_vector_noPBC(n_p),y_vector_noPBC(n_p),&
                                  z_vector_noPBC(n_p)                       ! componentes del vector posición sin PBC
    real(dp),    intent(inout) :: wxx_matrix(n_p,tau_max_corr),&            ! matrices de acumulación
                                  wyy_matrix(n_p,tau_max_corr),&
                                  wzz_matrix(n_p,tau_max_corr)
    real(dp),    intent(inout) :: sum_wxx_vector(tau_max_corr),&            ! vectores de sumas auxiliares
                                  sum_wyy_vector(tau_max_corr),&
                                  sum_wzz_vector(tau_max_corr)
    integer(sp),    intent(inout) :: counter_data(tau_max_corr)                ! contador de datos
    integer(sp), intent(inout) :: counter                                   ! contador de entradas
    
    integer(sp)            :: i,j
    integer(sp)            :: tau_corr_0,tau_corr_t ! tiempos de correlación
    ! NOTA: CONDICIÓN QUE SE DEBE CUMPLIR
    ! nmax_tau_corr_0 < (time_eq+time_run)/(ensamble_step*tau_max_corr)
    integer(sp), parameter :: nmax_tau_corr_0=5_sp ! maximo número de tau_corr_0 que almacenamos

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
            counter_data(j)=counter_data(j)+1_sp
        end do
    end if
end subroutine mean_squared_displacement