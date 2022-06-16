! make clean && make md_lj_canonical_ensamble_02.o && ./md_lj_canonical_ensamble_02.o
program md_lj_canonical_ensamble_02
    use module_precision;use module_md_lennard_jones
    implicit none
    integer(sp), parameter   :: n_p=500_sp                             ! cantidad de partículasa
    real(dp),    parameter   :: delta_time=0.005_dp                    ! paso temporal
    integer(sp), parameter   :: time_eq=2000_sp,&                      ! pasos de equilibración
                                time_run=1000_sp                       ! pasos de evolucion en el estado estacionario
    real(dp),    parameter   :: T_adim_ref=1.0_dp                      ! temperatura de referencia adimensional
    real(dp),    parameter   :: r_cutoff=2.5_dp,mass=1._dp             ! radio de corte de interacciones y masa     
    real(dp),    allocatable :: x_vector(:),y_vector(:),z_vector(:)    ! componentes de las posiciones/particula
    real(dp),    allocatable :: vx_vector(:),vy_vector(:),vz_vector(:) ! componentes de la velocidad/particula
    real(dp),    allocatable :: force_x(:),force_y(:),force_z(:)       ! componentes de la fuerza/particula
    integer(sp)              :: i,index,istat                                ! loop index
    real(dp)                 :: T_adim                                 ! Temperatura
    real(dp)                 :: vx_mc,vy_mc,vz_mc                      ! componentes de la velocidad del centro de masas
    real(dp)                 :: time,time_end,time_start               ! tiempos de CPU
    real(dp),    parameter   :: density=0.8_dp                         ! densidad (particulas/volumen)
    !real(dp),    parameter   :: density=1.2_dp                         ! densidad (particulas/volumen)

    ! DESCOMENTAR PARA density=0.8
    open(10,file='../results/structure_function_rho1.dat',status='replace',action='write',iostat=istat)
    ! DESCOMENTAR PARA density=1.2
    !open(10,file='../results/structure_function_rho2.dat',status='replace',action='write',iostat=istat)
    if (istat/=0) write(*,*) 'ERROR! istat(11file) = ',istat
    24 format(E12.4,x,E12.4);25 format(A12,x,A12)
    write(10,25) 'time','S(k,t)'

    call cpu_time(time_start)

    allocate(x_vector(n_p),y_vector(n_p),z_vector(n_p))
    x_vector(:)=0._dp;y_vector(:)=0._dp;z_vector(:)=0._dp

    ! generamos configuración inicial (FCC structure)
    call initial_lattice_configuration(n_p,density,x_vector,y_vector,z_vector,2)
    write(10,24) 0._dp,static_structure_factor(n_p,density,x_vector,y_vector,z_vector)

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
        vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force_x,force_y,force_z)
        ! velocity center of mass to zero
        vx_mc=sum(vx_vector(:))*(1._dp/real(n_p,dp));vx_vector(:)=(vx_vector(:)-vx_mc)
        vy_mc=sum(vy_vector(:))*(1._dp/real(n_p,dp));vy_vector(:)=(vy_vector(:)-vy_mc)
        vz_mc=sum(vz_vector(:))*(1._dp/real(n_p,dp));vz_vector(:)=(vz_vector(:)-vz_mc)
        T_adim=temperature(n_p,mass,vx_vector,vy_vector,vz_vector)
        time=real(i,dp)*delta_time
        write(10,24) time,static_structure_factor(n_p,density,x_vector,y_vector,z_vector)
    end do

    ! ESTACIONARIO
    time=0._dp
    do i=1,time_run
        index=index+1
        write(*,*) 'paso temporal =',index,' de',time_eq+time_run
        call rescaling_velocities(n_p,vx_vector,vy_vector,vz_vector,T_adim_ref,mass)
        call velocity_verlet(n_p,x_vector,y_vector,z_vector,&
        vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force_x,force_y,force_z)
        time=(real(time_eq,dp)+real(i,dp))*delta_time
        write(10,24) time,static_structure_factor(n_p,density,x_vector,y_vector,z_vector)
    end do
    close(10)

    deallocate(x_vector,y_vector,z_vector)
    deallocate(vx_vector,vy_vector,vz_vector)
    deallocate(force_x,force_y,force_z)

    call cpu_time(time_end)
    write(*,*) 'elapsed time = ',time_end-time_start,'[s]'
end program md_lj_canonical_ensamble_02