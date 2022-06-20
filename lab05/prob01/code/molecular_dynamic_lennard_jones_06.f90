! Inciso g)
! make clean && make molecular_dynamic_lennard_jones_06.o && ./molecular_dynamic_lennard_jones_06.o
program molecular_dynamic_lennard_jones_06
    use module_precision;use module_md_lennard_jones
    implicit none
    integer(sp), parameter   :: nc_max=8_sp                            ! factor máximo paa definir número de partículas
    integer(sp), parameter   :: time_eq=1000_sp,time_scal=50_sp,&      ! pasos de equilibración y de escaleo de veloc.
                                time_run=10_sp                         ! pasos de evolucion en el estado estacionario
    real(dp),    parameter   :: T_adim_ref=1.1_dp                      ! temperatura de referencia adimensional
    real(dp),    parameter   :: density=0.8_dp                         ! densidad (particulas/volumen)
    real(dp),    parameter   :: mass=1._dp                             ! masa
    real(dp),    parameter   :: delta_time=0.005_dp                    ! paso temporal     
    real(dp),    allocatable :: x_vector(:),y_vector(:),z_vector(:)    ! componentes de las posiciones/particula
    real(dp),    allocatable :: vx_vector(:),vy_vector(:),vz_vector(:) ! componentes de la velocidad/particula
    real(dp),    allocatable :: force_x(:),force_y(:),force_z(:)       ! componentes de la fuerza/particula
    integer(sp)              :: n_p                                    ! cantidad de partículas
    integer(sp)              :: i,j,istat,index                        ! loop index
    real(dp)                 :: r_cutoff,L                             ! radio de corte de interacciones
    real(dp)                 :: T_adim                                 ! temperatura adimensional
    real(dp)                 :: vx_mc,vy_mc,vz_mc                      ! componentes de la velocidad del centro de masas
    real(dp)                 :: time_end,time_start                    ! tiempos de CPU

    open(10,file='../results/num_particles_vs_cpu_time.dat',status='replace',action='write',iostat=istat)
    if (istat/=0) write(*,*) 'ERROR! istat(10file) = ',istat
    11 format(E12.4,x,I12);12 format(A12,x,A12)
    write(10,12) 'CPU_time','n_p'

    do j=1,nc_max
        call cpu_time(time_start)
        write(*,*) 'corrida=',j,' de ',nc_max
        n_p=4_sp*j*j*j
        L=(real(n_p,dp)*(1._dp/density))**(1._dp/3._dp)
        r_cutoff=L*0.3_dp ! definimos r_cutoff < L/2

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

        index=0
        ! TRANSITORIO
        do i=1,time_eq
            index=index+1
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
        do i=1,time_run
            index=index+1
            call velocity_verlet(n_p,x_vector,y_vector,z_vector,&
            vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force_x,force_y,force_z)
        end do
        write(*,*) 'termino el estacionario'

        deallocate(x_vector,y_vector,z_vector)
        deallocate(vx_vector,vy_vector,vz_vector)
        deallocate(force_x,force_y,force_z)
        call cpu_time(time_end)
        write(10,11) time_end-time_start,n_p
    end do
    close(10)
end program molecular_dynamic_lennard_jones_06