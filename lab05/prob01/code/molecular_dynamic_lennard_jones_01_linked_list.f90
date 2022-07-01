! make clean && make molecular_dynamic_lennard_jones_01_linked_list.o && ./molecular_dynamic_lennard_jones_01_linked_list.o
program molecular_dynamic_lennard_jones_01_linked_list
    use module_precision;use module_md_lennard_jones
    implicit none
    integer(sp), parameter   :: n_p=256_sp                             ! cantidad de partículasa
    real(dp),    parameter   :: delta_time=0.005_dp                    ! paso temporal
    integer(sp), parameter   :: time_eq=1000_sp,time_scal=50_sp,&      ! pasos de equilibración y de escaleo de veloc.
                                time_run=1000_sp                       ! pasos de evolucion en el estado estacionario
    real(dp),    parameter   :: T_adim_ref=1.1_dp                      ! temperatura de referencia adimensional
    real(dp),    parameter   :: density=0.8_dp                         ! densidad (particulas/volumen)
    real(dp),    parameter   :: r_cutoff=2.5_dp,mass=1._dp             ! radio de corte de interacciones y masa     
    real(dp),    allocatable :: x_vector(:),y_vector(:),z_vector(:),&  ! componentes de las posiciones/particula
                                x_vector_noPBC(:),y_vector_noPBC(:),&
                                z_vector_noPBC(:)
    real(dp),    allocatable :: vx_vector(:),vy_vector(:),vz_vector(:) ! componentes de la velocidad/particula
    real(dp),    allocatable :: force_x(:),force_y(:),force_z(:)       ! componentes de la fuerza/particula
    integer(sp)              :: i,istat                          ! loop index
    real(dp)                 :: U_adim,Ec_adim,time,press,v_mc,T_adim  ! observables
    real(dp)                 :: vx_mc,vy_mc,vz_mc                      ! componentes de la velocidad del centro de masas
    real(dp)                 :: time_end,time_start                    ! tiempos de CPU
    logical,     parameter   :: energies_switch=.true.                 ! escribir energías en el estado estacionario
    ! variable para aplicar linked list
    integer(sp), allocatable :: map(:),list(:),head(:)
    integer(sp)              :: m
    real(dp)                 :: L

    call cpu_time(time_start)
    22 format(5(E12.4,x),E12.4);23 format(5(A12,x),A12)

    allocate(x_vector(n_p),y_vector(n_p),z_vector(n_p))
    x_vector(:)=0._dp;y_vector(:)=0._dp;z_vector(:)=0._dp
    allocate(x_vector_noPBC(n_p),y_vector_noPBC(n_p),z_vector_noPBC(n_p))
    x_vector_noPBC(:)=0._dp;y_vector_noPBC(:)=0._dp;z_vector_noPBC(:)=0._dp

    L=(n_p*(1._dp/density))**(1._dp/3._dp)
    m=int((1._dp/int(L*(1._dp/r_cutoff),sp))*L,sp)

    allocate(map(13*m*m*m),list(n_p),head(m*m*m))
    map(:)=0;list(:)=0;head(:)=0
    call maps(m,map)

    ! generamos configuración inicial (FCC structure)
    call initial_lattice_configuration(n_p,density,x_vector,y_vector,z_vector,2)
    x_vector_noPBC(:)=x_vector(:);y_vector_noPBC(:)=y_vector(:);z_vector_noPBC(:)=z_vector(:)
    call links(n_p,m,L,head,list,x_vector,y_vector,z_vector)

    allocate(vx_vector(n_p),vy_vector(n_p),vz_vector(n_p))
    vx_vector(:)=0._dp;vy_vector(:)=0._dp;vz_vector(:)=0._dp
    call md_initial_parameters(n_p,x_vector,y_vector,z_vector,&
    vx_vector,vy_vector,vz_vector,T_adim_ref,delta_time,density,mass)

    ! computamos fuerzas en el tiempo inicial
    allocate(force_x(n_p),force_y(n_p),force_z(n_p))
    force_x(:)=0._dp;force_y(:)=0._dp;force_z(:)=0._dp
    ! call f_lj_total(x_vector,y_vector,z_vector,r_cutoff,n_p,density,&
    !     force_x,force_y,force_z)
    call f_lj_total_linkedlist(x_vector,y_vector,z_vector,r_cutoff,n_p,density,&
            force_x,force_y,force_z,m,map,list,head)

    ! TRANSITORIO
    do i=1,time_eq
        write(*,*) i
        if (mod(i,time_scal)==0_sp) call rescaling_velocities(n_p,vx_vector,vy_vector,vz_vector,T_adim_ref,mass)

        ! call velocity_verlet(n_p,x_vector,y_vector,z_vector,&
        ! x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
        ! vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force_x,force_y,force_z)

        call velocity_verlet_linked_list(n_p,x_vector,y_vector,z_vector,&
                x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
                vx_vector,vy_vector,vz_vector,&
                delta_time,mass,r_cutoff,density,force_x,force_y,force_z,&
                m,map,list,head)

        ! velocity center of mass to zero
        vx_mc=sum(vx_vector(:))*(1._dp/real(n_p,dp));vx_vector(:)=(vx_vector(:)-vx_mc)
        vy_mc=sum(vy_vector(:))*(1._dp/real(n_p,dp));vy_vector(:)=(vy_vector(:)-vy_mc)
        vz_mc=sum(vz_vector(:))*(1._dp/real(n_p,dp));vz_vector(:)=(vz_vector(:)-vz_mc)
        T_adim=temperature(n_p,mass,vx_vector,vy_vector,vz_vector)
    end do

    ! ESTACIONARIO
    time=0._dp

    ! ESCRIBIMOS DATOS
    if (energies_switch.eqv..true.) then
        open(12,file='../results/result_03_linked_list.dat',status='replace',action='write',iostat=istat)
        write(*,*) 'istat(12file) = ',istat;write(12,23) 'time','pot_ergy','kin_ergy','v_mc','press','T_adim'
    end if

    U_adim=u_lj_total(n_p,x_vector,y_vector,z_vector,r_cutoff,density)
    Ec_adim=kinetic_ergy_total(n_p,vx_vector,vy_vector,vz_vector,mass)
    press=pressure(n_p,density,mass,r_cutoff,x_vector,y_vector,z_vector,&
    vx_vector,vy_vector,vz_vector)
    T_adim=temperature(n_p,mass,vx_vector,vy_vector,vz_vector)
    v_mc=sqrt(vx_mc*vx_mc+vy_mc*vy_mc+vz_mc*vz_mc)

    ! ESCRIBIMOS DATOS
    if (energies_switch.eqv..true.) then
        write(12,22) time,U_adim*(1._dp/real(n_p,dp)),Ec_adim*(1._dp/real(n_p,dp)),v_mc,press,T_adim
    end if


    do i=1,time_run
        write(*,*) time_eq+i

        ! call velocity_verlet(n_p,x_vector,y_vector,z_vector,&
        ! x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
        ! vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force_x,force_y,force_z)

        call velocity_verlet_linked_list(n_p,x_vector,y_vector,z_vector,&
                x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
                vx_vector,vy_vector,vz_vector,&
                delta_time,mass,r_cutoff,density,force_x,force_y,force_z,&
                m,map,list,head)

        ! velocity center of mass to zero
        vx_mc=sum(vx_vector(:))*(1._dp/real(n_p,dp));vx_vector(:)=(vx_vector(:)-vx_mc)
        vy_mc=sum(vy_vector(:))*(1._dp/real(n_p,dp));vy_vector(:)=(vy_vector(:)-vy_mc)
        vz_mc=sum(vz_vector(:))*(1._dp/real(n_p,dp));vz_vector(:)=(vz_vector(:)-vz_mc)
        U_adim=u_lj_total(n_p,x_vector,y_vector,z_vector,r_cutoff,density)
        Ec_adim=kinetic_ergy_total(n_p,vx_vector,vy_vector,vz_vector,mass)
        press=pressure(n_p,density,mass,r_cutoff,x_vector,y_vector,z_vector,&
        vx_vector,vy_vector,vz_vector)
        T_adim=temperature(n_p,mass,vx_vector,vy_vector,vz_vector)
        v_mc=sqrt(vx_mc*vx_mc+vy_mc*vy_mc+vz_mc*vz_mc)
        time=real(i,dp)*delta_time

        ! ESCRIBIMOS DATOS
        if (energies_switch.eqv..true.) then
            write(12,22) time,U_adim*(1._dp/real(n_p,dp)),Ec_adim*(1._dp/real(n_p,dp)),v_mc,press,T_adim
        end if
    end do
    if (energies_switch.eqv..true.) close(12)

    deallocate(x_vector,y_vector,z_vector)
    deallocate(vx_vector,vy_vector,vz_vector)
    deallocate(force_x,force_y,force_z)
    deallocate(map,list,head)

    call cpu_time(time_end)
    write(*,*) 'elapsed time = ',time_end-time_start,'[s]'
end program molecular_dynamic_lennard_jones_01_linked_list