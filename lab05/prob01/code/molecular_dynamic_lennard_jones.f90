! make clean && make molecular_dynamic_lennard_jones.o && ./molecular_dynamic_lennard_jones.o
program molecular_dynamic_lennard_jones
    use module_precision;use module_md_lennard_jones
    implicit none
    integer(sp), parameter   :: n_p=125_sp                             ! sitios de red por dimension
    real(dp),    parameter   :: delta_time=0.005_dp
    integer(sp), parameter   :: time_eq=1000_sp,time_scal=50_sp,time_run=1000_sp
    real(dp),    parameter   :: T_adim_ref=1.1_dp                      ! temperatura de referencia adimensional
    real(dp),    parameter   :: density=0.8_dp
    real(dp),    parameter   :: r_cutoff=2.5_dp,mass=1._dp
    real(dp),    allocatable :: x_vector(:),y_vector(:),z_vector(:) ! vectores posicion
    real(dp),    allocatable :: vx_vector(:),vy_vector(:),vz_vector(:) ! vectores velocidad
    real(dp),    allocatable :: force(:)
    integer(sp)              :: i,istat
    real(dp)                 :: U_adim,Ec_adim,time,press,v_mc
    real(dp)                 :: vx_mc,vy_mc,vz_mc
    real(dp)                 :: time_end,time_start

    call cpu_time(time_start)
    20 format(2(E12.4,x),E12.4);21 format(2(A12,x),A12)
    22 format(4(E12.4,x),E12.4);23 format(4(A12,x),A12)

    allocate(x_vector(n_p),y_vector(n_p),z_vector(n_p))
    x_vector(:)=0._dp;y_vector(:)=0._dp;z_vector(:)=0._dp

    ! generamos configuración inicial (FCC structure)
    call initial_lattice_configuration(n_p,density,x_vector,y_vector,z_vector,1)

    ! ESCRIBIMOS DATOS
    open(10,file='../results/result_01.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat;write(10,21) 'rx_fcc','ry_fcc','rz_fcc'
    do i=1,n_p;write(10,20) x_vector(i),y_vector(i),z_vector(i);end do;close(10)

    allocate(vx_vector(n_p),vy_vector(n_p),vz_vector(n_p))
    vx_vector(:)=0._dp;vy_vector(:)=0._dp;vz_vector(:)=0._dp
    call md_initial_parameters(n_p,x_vector,y_vector,z_vector,&
    vx_vector,vy_vector,vz_vector,T_adim_ref,delta_time,density,mass)

    ! computamos fuerzas en el tiempo inicial
    allocate(force(n_p));force(:)=0._dp
    call f_lj_total(x_vector,y_vector,z_vector,r_cutoff,n_p,density,force)

    ! TRANSITORIO
    do i=1,time_eq
        write(*,*) i
        if (mod(i,time_scal)==0_sp) call rescaling_velocities(n_p,vx_vector,vy_vector,vz_vector,T_adim_ref,mass)
        call velocity_verlet(n_p,x_vector,y_vector,z_vector,&
        vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force)
        ! velocity center of mass to zero
        vx_mc=sum(vx_vector(:))*(1._dp/real(n_p,dp));vx_vector(:)=(vx_vector(:)-vx_mc)
        vy_mc=sum(vy_vector(:))*(1._dp/real(n_p,dp));vy_vector(:)=(vy_vector(:)-vy_mc)
        vz_mc=sum(vz_vector(:))*(1._dp/real(n_p,dp));vz_vector(:)=(vz_vector(:)-vz_mc)
    end do

    ! ESTACIONARIO
    time=0._dp
    open(12,file='../results/result_03.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(12file) = ',istat;write(12,23) 'time','pot_ergy','kin_ergy','v_mc','press'
    U_adim=u_lj_total(n_p,x_vector,y_vector,z_vector,r_cutoff,density)
    Ec_adim=kinetic_ergy_total(n_p,vx_vector,vy_vector,vz_vector,mass)
    press=pressure(n_p,density,mass,r_cutoff,x_vector,y_vector,z_vector,&
    vx_vector,vy_vector,vz_vector)
    v_mc=sqrt(vx_mc*vx_mc+vy_mc*vy_mc+vz_mc*vz_mc)
    write(12,22) time,U_adim*(1._dp/real(n_p,dp)),Ec_adim*(1._dp/real(n_p,dp)),v_mc,press
    do i=1,time_run
        write(*,*) time_eq+i
        call velocity_verlet(n_p,x_vector,y_vector,z_vector,&
        vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force)
        ! velocity center of mass to zero
        vx_mc=sum(vx_vector(:))*(1._dp/real(n_p,dp));vx_vector(:)=(vx_vector(:)-vx_mc)
        vy_mc=sum(vy_vector(:))*(1._dp/real(n_p,dp));vy_vector(:)=(vy_vector(:)-vy_mc)
        vz_mc=sum(vz_vector(:))*(1._dp/real(n_p,dp));vz_vector(:)=(vz_vector(:)-vz_mc)
        U_adim=u_lj_total(n_p,x_vector,y_vector,z_vector,r_cutoff,density)
        Ec_adim=kinetic_ergy_total(n_p,vx_vector,vy_vector,vz_vector,mass)
        press=pressure(n_p,density,mass,r_cutoff,x_vector,y_vector,z_vector,&
        vx_vector,vy_vector,vz_vector)
        v_mc=sqrt(vx_mc*vx_mc+vy_mc*vy_mc+vz_mc*vz_mc)
        time=real(i,dp)*delta_time
        write(12,22) time,U_adim*(1._dp/real(n_p,dp)),Ec_adim*(1._dp/real(n_p,dp)),v_mc,press
    end do
    close(12)

    ! ESCRIBIMOS DATOS
    open(11,file='../results/result_02.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(11file) = ',istat;write(11,21) 'rx_fcc','ry_fcc','rz_fcc'
    do i=1,n_p;write(11,20) x_vector(i),y_vector(i),z_vector(i);end do;close(11)

    deallocate(x_vector,y_vector,z_vector)
    deallocate(vx_vector,vy_vector,vz_vector)
    deallocate(force)
    call cpu_time(time_end)
    write(*,*) 'elapsed time = ',time_end-time_start,'[s]'
end program molecular_dynamic_lennard_jones