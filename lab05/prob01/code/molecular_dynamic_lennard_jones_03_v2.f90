! make clean && make molecular_dynamic_lennard_jones_03_v2.o && ./molecular_dynamic_lennard_jones_03_v2.o
program molecular_dynamic_lennard_jones_03_v2
    use module_precision;use module_md_lennard_jones
    implicit none
    integer(sp), parameter   :: n_p=256_sp                             ! cantidad de partículasa
    real(dp),    parameter   :: T_adim_ref=1.1_dp                      ! temperatura de referencia adimensional
    real(dp),    parameter   :: density=0.8_dp                         ! densidad (particulas/volumen)
    real(dp),    parameter   :: r_cutoff=2.5_dp,mass=1._dp             ! radio de corte de interacciones y masa     
    real(dp),    allocatable :: x_vector(:),y_vector(:),z_vector(:)    ! componentes de las posiciones/particula
    real(dp),    allocatable :: vx_vector(:),vy_vector(:),vz_vector(:) ! componentes de la velocidad/particula
    real(dp),    allocatable :: force_x(:),force_y(:),force_z(:)       ! componentes de la fuerza/particula
    integer(sp)              :: time_eq,time_scal,&                    ! pasos de equilibración y de escaleo de veloc.
                                time_run                               ! pasos de evolucion en el estado estacionario
    real(dp)                 :: delta_time                             ! paso temporal
    integer(sp)              :: i,j,istat,index                        ! loop index
    real(dp)                 :: U_adim,U_med,var_U,err_U
    real(dp)                 :: Ec_adim,Ec_med,var_Ec,err_Ec
    real(dp)                 :: Etot_adim,Etot_med,var_Etot,err_Etot
    real(dp)                 :: T_adim
    real(dp)                 :: s1_U,s2_U
    real(dp)                 :: s1_Ec,s2_Ec
    real(dp)                 :: s1_Etot,s2_Etot
    real(dp)                 :: vx_mc,vy_mc,vz_mc                      ! componentes de la velocidad del centro de masas
    real(dp)                 :: time_end,time_start                    ! tiempos de CPU

    call cpu_time(time_start)

    allocate(x_vector(n_p),y_vector(n_p),z_vector(n_p))
    x_vector(:)=0._dp;y_vector(:)=0._dp;z_vector(:)=0._dp

    ! generamos configuración inicial (FCC structure)
    call initial_lattice_configuration(n_p,density,x_vector,y_vector,z_vector,2)

    allocate(vx_vector(n_p),vy_vector(n_p),vz_vector(n_p))
    vx_vector(:)=0._dp;vy_vector(:)=0._dp;vz_vector(:)=0._dp
    call md_initial_parameters(n_p,x_vector,y_vector,z_vector,&
    vx_vector,vy_vector,vz_vector,T_adim_ref,delta_time,density,mass)

    ! computamos fuerzas en el tiempo inicial
    allocate(force_x(n_p),force_y(n_p),force_z(n_p))
    force_x(:)=0._dp;force_y(:)=0._dp;force_z(:)=0._dp
    call f_lj_total(x_vector,y_vector,z_vector,r_cutoff,n_p,density,force_x,force_y,force_z)
 

    open(10,file='../results/enrgy_fluct_relations.dat',status='replace',action='write',iostat=istat)
    if (istat/=0) write(*,*) 'ERROR! istat(10file) = ',istat
    11 format(E12.4,x,E12.4);12 format(A12,x,A12)
    write(10,12) 'delta_time','var_Etot/var_Ec'

    delta_time=0.02_dp
    do j=1,5
        ! definimos los pasos temporales y pasos de MD
        time_eq=5._dp*(1._dp/delta_time)
        time_scal=0.25_dp*(1._dp/delta_time)
        time_run=5._dp*(1._dp/delta_time)

        index=0
        ! TRANSITORIO
        do i=1,time_eq
            index=index+1
            write(*,*) 'paso temporal =',index,' de',time_eq+time_run
            if (mod(i,time_scal)==0_sp) call rescaling_velocities(n_p,vx_vector,vy_vector,vz_vector,T_adim_ref,mass)
            call velocity_verlet(n_p,x_vector,y_vector,z_vector,&
            vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force_x,force_y,force_z)
            ! velocity center of mass to zero
            vx_mc=sum(vx_vector(:))*(1._dp/real(n_p,dp));vx_vector(:)=(vx_vector(:)-vx_mc)
            vy_mc=sum(vy_vector(:))*(1._dp/real(n_p,dp));vy_vector(:)=(vy_vector(:)-vy_mc)
            vz_mc=sum(vz_vector(:))*(1._dp/real(n_p,dp));vz_vector(:)=(vz_vector(:)-vz_mc)
            T_adim=temperature(n_p,mass,vx_vector,vy_vector,vz_vector)
        end do

        write(*,*) 'termino el transitorio'

        ! ESTACIONARIO
        U_med=0._dp
        Ec_med=0._dp
        Etot_med=0._dp

        s1_U=0._dp;s2_U=0._dp
        s1_Ec=0._dp;s2_Ec=0._dp
        s1_Etot=0._dp;s2_Etot=0._dp
        
        do i=1,time_run
            index=index+1
            write(*,*) 'paso temporal =',index,' de',time_eq+time_run
            call velocity_verlet(n_p,x_vector,y_vector,z_vector,&
            vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force_x,force_y,force_z)

            U_adim=u_lj_total(n_p,x_vector,y_vector,z_vector,r_cutoff,density)
            Ec_adim=kinetic_ergy_total(n_p,vx_vector,vy_vector,vz_vector,mass)
            Etot_adim=(U_adim+Ec_adim)*(1._dp/real(n_p,dp))

            ! computamos 1er y 2do momento
            s1_U=s1_U+U_adim*(1._dp/real(n_p,dp));s2_U=s2_U+U_adim*U_adim*(1._dp/real(n_p*n_p,dp))
            s1_Ec=s1_Ec+Ec_adim*(1._dp/real(n_p,dp));s2_Ec=s2_Ec+Ec_adim*Ec_adim*(1._dp/real(n_p*n_p,dp))
            s1_Etot=s1_Etot+Etot_adim;s2_Etot=s2_Etot+Etot_adim*Etot_adim

            ! computamos valores medios (mejor a mayor paso evolucionado)
            U_med=s1_U*(1._dp/real(i,dp))
            Ec_med=s1_Ec*(1._dp/real(i,dp))
            Etot_med=s1_Etot*(1._dp/real(i,dp))

            ! computamos varianzas (mejor a mayor paso evolucionado)
            var_U=(real(i,dp)*s2_U-s1_U*s1_U)*(1._dp/real(i*i,dp))
            var_Ec=(real(i,dp)*s2_Ec-s1_Ec*s1_Ec)*(1._dp/real(i*i,dp))
            var_Etot=(real(i,dp)*s2_Etot-s1_Etot*s1_Etot)*(1._dp/real(i*i,dp))

        end do

        write(*,*) 'termino el estacionario'

        ! computamos errores en el último paso
        err_U=(var_U*0.25_dp)*(1._dp/real(time_eq-1,dp))
        err_Ec=(var_Ec*0.25_dp)*(1._dp/real(time_eq-1,dp))
        err_Etot=(var_Etot*0.25_dp)*(1._dp/real(time_eq-1,dp))

        write(10,11) delta_time,var_Etot*(1._dp/var_Ec)

        delta_time=delta_time*0.5_dp
    end do

    close(10)

    deallocate(x_vector,y_vector,z_vector)
    deallocate(vx_vector,vy_vector,vz_vector)
    deallocate(force_x,force_y,force_z)

    call cpu_time(time_end)
    write(*,*) 'elapsed time = ',time_end-time_start,'[s]'
end program molecular_dynamic_lennard_jones_03_v2