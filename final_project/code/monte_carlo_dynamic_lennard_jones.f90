! make clean && make monte_carlo_dynamic_lennard_jones.o && ./monte_carlo_dynamic_lennard_jones.o
program monte_carlo_dynamic_lennard_jones
    use module_precision;use module_mc_lennard_jones
    implicit none
    ! VARIABLES GENERALES
    integer(sp), parameter   :: n_p=256_sp                                  ! cantidad de partículasa
    real(dp),    parameter   :: delta_time=0.05_dp
    integer(sp), parameter   :: MC_step_eq=0_sp!1_sp!                       ! monte carlo step para equilibración (transitorio)
    integer(sp), parameter   :: MC_step_run=100_sp!1000_sp!0_sp!           ! monete carlo step para corrida (estacionario)
    real(dp),    parameter   :: T_adim_ref=0.75_dp                          ! temperatura de referencia adimensional
    real(dp),    parameter   :: r_cutoff=2.5_dp,mass=1._dp                  ! radio de corte de interacciones y masa     
    real(dp)                 :: delta_x,delta_y,delta_z
    real(dp),    allocatable :: x_vector(:),y_vector(:),z_vector(:)         ! componentes de las posiciones/particula
    integer(sp)              :: i,j,k,index,istat                           ! loop index
    real(dp)                 :: time,L                                        ! tiempo de evolución
    ! VARIABLES PARA COMPUTAR FUERZAS CON LINKED LIST
    integer(sp), parameter   :: m=3_sp                                      ! numero de celdas por dimensión 
    integer(sp), allocatable :: map(:),list(:),head(:)
    integer(sp), parameter   :: linkedlist_type=1_sp                    ! simular con(1)/sin(0) linkedlist
    ! VARIABLES PARA REALIZAR TERMALIZACIÓN GRADUAL (EVITANDO QUENCHIN)
    real(dp),   parameter    :: T_adim_min=T_adim_ref*0.1_dp,&              ! Temperatura mínima
                                T_adim_max=T_adim_ref*0.95_dp               ! Temperatura máxima
    integer(sp), parameter   :: n_Temp=10_sp                                ! cantidad de temperaturas simuladas
    real(dp),    parameter   :: step_Temp=abs(T_adim_max-T_adim_min)*&      ! paso de variación de temperaturas
                                             (1._dp/real(n_Temp,dp))
    real(dp)                 :: Tadim
    ! VARIABLES PARA COMPUTAR TIEMPO TRANSCURRIDO DE CPU
    real(dp)                 :: time_end,time_start                          ! tiempos de CPU
    ! VARIABLES PARA ACTIVAR/DESACTIVAR ESCRITURA DE DATOS
    logical                  :: pressure_switch=.false.,&                    ! presión vs densidad
                                structure_factor_switch=.false.,&            ! factor de estructura vs densidad
                                diffusion_coeff_switch=.false.,&             ! coeficiente de difusión vs densidad
                                energie_switch=.true.                        ! energía interna
    ! VARIABLES PARA REALIZAR BARRIDO DE DENSIDADES
    real(dp),    parameter   :: density_min=0.8_dp,density_max=0.8_dp       ! rango de densidades
    integer(sp), parameter   :: n_density=1_sp                              ! cantidad de densidades simuladas
    real(dp),    parameter   :: step_density=abs(density_max-density_min)*&
                                (1._dp/real(n_density,dp))                  ! paso de variación de densidades
    real(dp)                 :: density                                     ! densidad (particulas/volumen)
    ! VARIABLES PARA COMPUTAR ENERGÍA INTERNA
    real(dp)                 :: Uadim,Uadim_med,var_Uadim,err_Uadim         ! energía interna adimensional
    real(dp)                 :: s1_Uadim,s2_Uadim                           ! 1er y 2do momento
    ! VARIABLES PARA COMPUTAR PRESIÓN OSMÓTICA
    real(dp)                 :: press,press_med,var_press,err_press         ! presión osmótica
    real(dp)                 :: s1_press,s2_press                           ! 1er y 2do momento
    ! VARIABLES PARA COMPUTAR FACTOR DE ESTRUCTURA ESTÁTICO
    real(dp)                 :: Sk,Sk_med,var_Sk,err_Sk                     ! factor de estructura
    real(dp)                 :: s1_Sk,s2_Sk                                 ! 1er y 2do momento
    ! VARIABLES PARA COMPUTAR MSD Y COEFICIENTE DE DIFUSIÓN
    integer(sp), parameter   :: tau_max_corr=100_sp                         ! pasos maximos de correlación
    real(dp),    allocatable :: x_vector_noPBC(:),y_vector_noPBC(:),&       ! componentes de la posición sin PBC
                                z_vector_noPBC(:)
    real(dp),    allocatable :: wxx_matrix(:,:),wyy_matrix(:,:),&           ! matrices auxiliares para cálculo de msd
                                wzz_matrix(:,:)
    real(dp),    allocatable :: sum_wxx_vector(:),sum_wyy_vector(:),&       ! vectores auxiliares para cálculo de msd
                                sum_wzz_vector(:)
    integer(sp), allocatable :: counter_data(:)
    integer(sp)              :: counter
    real(dp)                 :: D,D_med,var_D,err_D
    real(dp)                 :: s1_D,s2_D
    real(dp)                 :: msd,msd_med,var_msd,err_msd
    real(dp)                 :: s1_msd,s2_msd

    ! Mensaje del progreso de la simulación
    print *, 'STARTING SIMULATION'

    ! comenzamos a medir tiempo de cpu
    call cpu_time(time_start)

    ! FORMATOS DE ESCRITURA A UTILIZAR
    20 format(2(E12.4,x),x,E12.4);21 format(I12,x,E12.4);22 format(4(E12.4,x),x,E12.4)
    ! APERTURA DE ARCHIVOS DE DATOS
    if (pressure_switch.eqv..true.) then
        open(10,file='../results/mcd_pressure_vs_density_T0.75.dat',status='replace',action='write',iostat=istat)
        if (istat/=0) write(*,*) 'ERROR! istat(10file) = ',istat
        write(10,'(2(A12,x),x,A12)') 'density','pressure','error';end if
    if (structure_factor_switch.eqv..true.) then
        open(11,file='../results/mcd_struct_factor_vs_density_T0.75.dat',status='replace',action='write',iostat=istat)
        if (istat/=0) write(*,*) 'ERROR! istat(11file) = ',istat
        write(11,'(2(A12,x),x,A12)') 'density','S(k)_med','error';end if
    if (diffusion_coeff_switch.eqv..true.) then
        open(12,file='../results/mcd_diffusion_vs_density_T0.75.dat',status='replace',action='write',iostat=istat)
        if (istat/=0) write(*,*) 'ERROR! istat(12file) = ',istat
        write(12,'(4(A12,x),x,A12)') 'density','D','error','msd','error';end if


    ! allocación de memoria
    allocate(x_vector(n_p),y_vector(n_p),z_vector(n_p))

    allocate(map(13*m*m*m),list(n_p),head(m*m*m))

    allocate(wxx_matrix(n_p,tau_max_corr),wyy_matrix(n_p,tau_max_corr),wzz_matrix(n_p,tau_max_corr))
    allocate(sum_wxx_vector(tau_max_corr),sum_wyy_vector(tau_max_corr),sum_wzz_vector(tau_max_corr))
    allocate(x_vector_noPBC(n_p),y_vector_noPBC(n_p),z_vector_noPBC(n_p))
    allocate(counter_data(tau_max_corr))
    

    ! barrido de densidades
    do j=1,n_density
        ! seteo de variables y parámetros
        x_vector(:)=0._dp;y_vector(:)=0._dp;z_vector(:)=0._dp

        map(:)=0_sp;list(:)=0_sp;head(:)=0_sp

        sum_wxx_vector(:)=0._dp;sum_wxx_vector(:)=0._dp;sum_wxx_vector(:)=0._dp
        x_vector_noPBC(:)=0._dp;y_vector_noPBC(:)=0._dp;z_vector_noPBC(:)=0._dp
        counter=0;counter_data(:)=0_sp

        ! definimos densidad en el rango [density_min;density_max]
        density=density_min+step_density*real(j,dp)

        ! generamos configuración inicial (FCC structure)
        call initial_lattice_configuration(n_p,density,x_vector,y_vector,z_vector,2)
        x_vector_noPBC(:)=x_vector(:);y_vector_noPBC(:)=y_vector(:);z_vector_noPBC(:)=z_vector(:)

        select case (linkedlist_type)
            case(1) ! simulation whit linked-list
                ! INICIALIZAMOS LISTA DE VECINOS
                L=(real(n_p,dp)*(1._dp/density))**(1._dp/3._dp)
                ! inicializo map (para usar linked list)
                call maps(m,map)
                ! INICIALIZAMOS LISTA DE VECINOS
                call links(n_p,m,L,head,list,x_vector,y_vector,z_vector)
                ! computamos energía interna en el tiempo inicial
                Uadim=u_lj_total_linkedlist(n_p,x_vector,y_vector,z_vector,&
                    r_cutoff,density,m,map,list,head)
            case(0) ! simulation whitout linked-list
                ! computamos energía interna en el tiempo inicial
                Uadim=u_lj_total(n_p,x_vector,y_vector,z_vector,r_cutoff,density)
        end select

        index=0_sp
        if (energie_switch.eqv..true.) then
            ! recordar medir energía vs tiempo para un único valor de densidad (CONTROL DE ESTABILIDAD)
            open(13,file='../results/mcd_energies.dat',status='replace',action='write',iostat=istat)
            if (istat/=0) write(*,*) 'ERROR! istat(13ile) = ',istat
            write(13,'(A12,x,A12)') 'MC_step','Uadim'
            write(13,21) index,Uadim
        end if

        ! computamos desplazamiendo cuadrático medio
        if (diffusion_coeff_switch.eqv..true.) then
            call mean_squared_displacement(n_p,x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,tau_max_corr,&
                wxx_matrix,wyy_matrix,wzz_matrix,sum_wxx_vector,sum_wyy_vector,sum_wzz_vector,&
                counter_data,counter)
        end if

        ! Termalizados hasta una temperatura anterior a la buscada
        do i=1,n_Temp
            Tadim=T_adim_min+step_Temp*real(i,dp)
            print*,i

            select case (linkedlist_type)
                case(1) ! simulation whit linked-list
                    ! calcular desplazamiento optimizados para acceptancia del 50%
                    call max_displacement_adjusting_linkedlist(n_p,x_vector,y_vector,z_vector,&
                        Uadim,Tadim,r_cutoff,density,delta_x,delta_y,delta_z,m,map,list,head)
                    call evolution_monte_carlo_linkedlist(n_p,x_vector,y_vector,z_vector,&
                        x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
                        Uadim,Tadim,r_cutoff,density,delta_x,delta_y,delta_z,m,map,list,head)
                case(0) ! simulation whitout linked-list
                    ! calcular desplazamiento optimizados para acceptancia del 50%
                    call max_displacement_adjusting(n_p,x_vector,y_vector,z_vector,&
                        Tadim,r_cutoff,density,delta_x,delta_y,delta_z)
                    call evolution_monte_carlo(n_p,x_vector,y_vector,z_vector,&
                        x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
                        Uadim,Tadim,r_cutoff,density,delta_x,delta_y,delta_z)
            end select
        end do

        select case (linkedlist_type)
            case(1) ! simulation whit linked-list
                ! calcular desplazamiento optimizados para acceptancia del 50%
                call max_displacement_adjusting_linkedlist(n_p,x_vector,y_vector,z_vector,&
                        Uadim,T_adim_ref,r_cutoff,density,delta_x,delta_y,delta_z,m,map,list,head)
            case(0) ! simulation whitout linked-list
                ! calcular desplazamiento optimizados para acceptancia del 50%
                call max_displacement_adjusting(n_p,x_vector,y_vector,z_vector,&
                        T_adim_ref,r_cutoff,density,delta_x,delta_y,delta_z)
        end select

        ! RÉGIMEN TRANSITORIO
        do i=1,MC_step_eq
            index=index+1
            ! Mensaje del progreso de la simulación
            print *, 'RUNNING...',(real(index,dp)/real((MC_step_eq+MC_step_run)*j,dp))*100._dp,'%',j

            select case (linkedlist_type)
                case(1) ! simulation whit linked-list
                    call evolution_monte_carlo_linkedlist(n_p,x_vector,y_vector,z_vector,&
                        x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
                        Uadim,Tadim,r_cutoff,density,delta_x,delta_y,delta_z,m,map,list,head)
                case(0) ! simulation whitout linked-list
                    call evolution_monte_carlo(n_p,x_vector,y_vector,z_vector,&
                        x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
                        Uadim,Tadim,r_cutoff,density,delta_x,delta_y,delta_z)
            end select
        end do

        ! RÉGIMEN ESTACIONARIO
        Uadim_med=0._dp;s1_Uadim=0._dp;s2_Uadim=0._dp
        press_med=0._dp;s1_press=0._dp;s2_press=0._dp
        Sk_med=0._dp;s1_Sk=0._dp;s2_Sk=0._dp

        do i=1,MC_step_run
            index=index+1
            ! Mensaje del progreso de la simulación
            print *, 'RUNNING...',(real(index,dp)/real((MC_step_eq+MC_step_run)*j,dp))*100._dp,'%',j

            select case (linkedlist_type)
                case(1) ! simulation whit linked-list
                    call evolution_monte_carlo_linkedlist(n_p,x_vector,y_vector,z_vector,&
                        x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
                        Uadim,Tadim,r_cutoff,density,delta_x,delta_y,delta_z,m,map,list,head)
                case(0) ! simulation whitout linked-list
                    call evolution_monte_carlo(n_p,x_vector,y_vector,z_vector,&
                        x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
                        Uadim,Tadim,r_cutoff,density,delta_x,delta_y,delta_z)
            end select

            if (pressure_switch.eqv..true.) then
                press=osmotic_pressure(n_p,density,mass,r_cutoff,x_vector,y_vector,z_vector)
                s1_press=s1_press+press;s2_press=s2_press+press*press
                press_med=s1_press*(1._dp/real(i,dp))
                var_press=(real(i,dp)*s2_press-s1_press*s1_press)*(1._dp/real(i*i,dp))
            end if
            if (structure_factor_switch.eqv..true.) then
                Sk=static_structure_factor(n_p,density,x_vector,y_vector,z_vector)
                s1_Sk=s1_Sk+Sk;s2_Sk=s2_Sk+Sk*Sk
                Sk_med=s1_Sk*(1._dp/real(i,dp))
                var_Sk=(real(i,dp)*s2_Sk-s1_Sk*s1_Sk)*(1._dp/real(i*i,dp))
            end if
            if (diffusion_coeff_switch.eqv..true.) then
                call mean_squared_displacement(n_p,x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,tau_max_corr,&
                    wxx_matrix,wyy_matrix,wzz_matrix,sum_wxx_vector,sum_wyy_vector,sum_wzz_vector,&
                    counter_data,counter)
            end if
            if (energie_switch.eqv..true.) then
                s1_Uadim=s1_Uadim+Uadim;s2_Uadim=s2_Uadim+Uadim*Uadim
                Uadim_med=s1_Uadim*(1._dp/real(i,dp))
                var_Uadim=(real(i,dp)*s2_Uadim-s1_Uadim*s1_Uadim)*(1._dp/real(i*i,dp))
                print*,Uadim_med
                write(13,21) index,Uadim_med
            end if
        end do

        if (energie_switch.eqv..true.) then
            ! computamos errores en el último paso
            err_Uadim=sqrt(var_Uadim*(1._dp/real(i-1,dp)))
            print*, 'Uadim_med=',Uadim_med,'+-',err_Uadim
        end if
        if (pressure_switch.eqv..true.) then
            ! computamos errores en el último paso
            err_press=sqrt(var_press*(1._dp/real(i-1,dp)))
            write(10,20) density,press_med,err_press
        end if
        if (structure_factor_switch.eqv..true.) then
            ! computamos errores en el último paso
            err_Sk=var_Sk*(1._dp/real(i-1,dp))
            write(11,20) density,Sk_med,err_Sk
        end if
        if (diffusion_coeff_switch.eqv..true.) then
            D_med=0._dp;s1_D=0._dp;s2_D=0._dp
            msd_med=0._dp;s1_msd=0._dp;s2_msd=0._dp
            do k=1,tau_max_corr
                if (counter_data(k)==0_sp) stop
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
            err_msd=sqrt(var_msd*(1._dp/real(tau_max_corr-1,dp)))
            err_D=sqrt(var_D*(1._dp/real(tau_max_corr-1,dp)))
            write(12,22) density,D_med,err_D,msd_med,err_msd
        end if
    end do

    if (pressure_switch.eqv..true.) close(10)
    if (structure_factor_switch.eqv..true.) close(11)
    if (diffusion_coeff_switch.eqv..true.) close(12)
    if (energie_switch.eqv..true.) close(13)

    ! liberamos memoria
    deallocate(x_vector,y_vector,z_vector)

    deallocate(map,list,head)

    deallocate(wxx_matrix,wyy_matrix,wzz_matrix)
    deallocate(sum_wxx_vector,sum_wyy_vector,sum_wzz_vector)
    deallocate(x_vector_noPBC,y_vector_noPBC,z_vector_noPBC)
    deallocate(counter_data)


    ! terminamos de medir tiempo de cpu
    call cpu_time(time_end)

    ! escribimos información de la simulación
    open(50,file='../results/mcd_data_run.dat',status='replace',action='write',iostat=istat)
        if (istat/=0) write(*,*) 'ERROR! istat(50file) = ',istat
        write(50,'(7(A12,x),x,A12)') 'cpu_time','delta_t','r_cutoff','T_ref',&
                                     'n_p','MC_step_eq','MC_step_run','tau_max_corr'
        write(50,'(4(E12.4,x),3(I12,x),I12)') time_end-time_start,delta_time,r_cutoff,T_adim_ref,&
                                              n_p,MC_step_eq,MC_step_run,tau_max_corr
    close(50)

    ! Mensaje del progreso de la simulación
    print *, 'FINISHING SIMULATION'  
end program monte_carlo_dynamic_lennard_jones

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
    integer(sp),    intent(inout) :: counter_data(tau_max_corr)                ! contador de datos
    integer(sp), intent(inout) :: counter                                   ! contador de entradas
    
    integer(sp)            :: i,j
    integer(sp)            :: tau_corr_0,tau_corr_t ! tiempos de correlación
    ! NOTA: CONDICIÓN QUE SE DEBE CUMPLIR
    ! nmax_tau_corr_0 < (time_eq+time_run)/tau_max_corr
    integer(sp), parameter :: nmax_tau_corr_0=5_sp ! maximo número de tau_corr_0 que almacenamos

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
            counter_data(j)=counter_data(j)+1_sp
        end do
    end if
end subroutine mean_squared_displacement