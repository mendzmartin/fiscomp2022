! make clean && make molecular_dynamic_lennard_jones.o && ./molecular_dynamic_lennard_jones.o
program molecular_dynamic_lennard_jones
    use module_precision;use module_md_lennard_jones
    implicit none
    integer(sp), parameter   :: n_p=256_sp                             ! cantidad de partículasa
    real(dp),    parameter   :: delta_time=0.005_dp                    ! paso temporal
    integer(sp), parameter   :: time_eq=1000_sp,time_scal=50_sp,&      ! pasos de equilibración y de escaleo de veloc.
                                time_run=1000_sp                       ! pasos de evolucion en el estado estacionario
    real(dp),    parameter   :: T_adim_ref=1.1_dp                      ! temperatura de referencia adimensional
    real(dp),    parameter   :: density=0.8_dp                         ! densidad (particulas/volumen)
    real(dp),    parameter   :: r_cutoff=2.5_dp,mass=1._dp             ! radio de corte de interacciones y masa     
    real(dp),    allocatable :: x_vector(:),y_vector(:),z_vector(:)    ! componentes de las posiciones/particula
    real(dp),    allocatable :: vx_vector(:),vy_vector(:),vz_vector(:) ! componentes de la velocidad/particula
    real(dp),    allocatable :: force_x(:),force_y(:),force_z(:)       ! componentes de la fuerza/particula
    integer(sp)              :: i,istat,index                          ! loop index
    real(dp)                 :: U_adim,Ec_adim,time,press,v_mc,T_adim  ! observables
    real(dp)                 :: vx_mc,vy_mc,vz_mc                      ! componentes de la velocidad del centro de masas
    real(dp)                 :: time_end,time_start                    ! tiempos de CPU
    logical                  :: movie_switch,fcc_init_switch,&
                                Tadim_trans_switch,energies_switch

    movie_switch       =.false. ! escribir peliculas en la caja
    fcc_init_switch    =.false. ! escribir estructura fcc inicial
    Tadim_trans_switch =.false. ! escribir tempereratura en el estado transitorio
    energies_switch    =.false.  ! escribir energías en el estado estacionario

    call cpu_time(time_start)
    22 format(5(E12.4,x),E12.4);23 format(5(A12,x),A12)

    allocate(x_vector(n_p),y_vector(n_p),z_vector(n_p))
    x_vector(:)=0._dp;y_vector(:)=0._dp;z_vector(:)=0._dp

    ! generamos configuración inicial (FCC structure)
    call initial_lattice_configuration(n_p,density,x_vector,y_vector,z_vector,2)

    ! ESCRIBIMOS DATOS
    if (fcc_init_switch.eqv..true.) then
        open(90,file='../results/fcc.dat',status='replace',action='write',iostat=istat)
        if (istat/=0) write(*,*) 'ERROR! istat(90file) = ',istat
        write(90,"(2(A12,x),A12)") 'rx_fcc','ry_fcc','rz_fcc'
        do i=1,n_p;write(90,"(2(E12.4,x),E12.4)") x_vector(i),y_vector(i),z_vector(i);end do;close(90)
    else if (movie_switch.eqv..true.) then
        index=10;call create_movie(index,x_vector,y_vector,z_vector,n_p)
    end if

    allocate(vx_vector(n_p),vy_vector(n_p),vz_vector(n_p))
    vx_vector(:)=0._dp;vy_vector(:)=0._dp;vz_vector(:)=0._dp
    call md_initial_parameters(n_p,x_vector,y_vector,z_vector,&
    vx_vector,vy_vector,vz_vector,T_adim_ref,delta_time,density,mass)

    ! computamos fuerzas en el tiempo inicial
    allocate(force_x(n_p),force_y(n_p),force_z(n_p))
    force_x(:)=0._dp;force_y(:)=0._dp;force_z(:)=0._dp
    call f_lj_total(x_vector,y_vector,z_vector,r_cutoff,n_p,density,force_x,force_y,force_z)

    ! ESCRIBIMOS DATOS
    if (fcc_init_switch.eqv..true.) then
        open(90,file='../results/init_force.dat',status='replace',action='write',iostat=istat)
        if (istat/=0) write(*,*) 'ERROR! istat(90file) = ',istat
        write(90,"(2(A12,x),A12)") 'fx','fy','fz'
        do i=1,n_p;write(90,"(2(E12.4,x),E12.4)") force_x(i),force_y(i),force_z(i);end do;close(90)
    else if (Tadim_trans_switch.eqv..true.) then
        open(90,file='../results/Tadim_transitorio.dat',status='replace',action='write',iostat=istat)
        if (istat/=0) write(*,*) 'ERROR! istat(90file) = ',istat
        write(90,"(A12,x,A12)") 'time','Tadim'
    end if

    ! TRANSITORIO
    do i=1,time_eq
        write(*,*) i
        if (mod(i,time_scal)==0_sp) call rescaling_velocities(n_p,vx_vector,vy_vector,vz_vector,T_adim_ref,mass)
        call velocity_verlet(n_p,x_vector,y_vector,z_vector,&
        vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force_x,force_y,force_z)
        ! velocity center of mass to zero
        vx_mc=sum(vx_vector(:))*(1._dp/real(n_p,dp));vx_vector(:)=(vx_vector(:)-vx_mc)
        vy_mc=sum(vy_vector(:))*(1._dp/real(n_p,dp));vy_vector(:)=(vy_vector(:)-vy_mc)
        vz_mc=sum(vz_vector(:))*(1._dp/real(n_p,dp));vz_vector(:)=(vz_vector(:)-vz_mc)
        T_adim=temperature(n_p,mass,vx_vector,vy_vector,vz_vector)
        if (Tadim_trans_switch.eqv..true.) write(90,"(E12.4,x,E12.4)") real(i,dp)*delta_time,T_adim
    end do
    if (Tadim_trans_switch.eqv..true.) close(90)

    ! ESTACIONARIO
    time=0._dp

    ! ESCRIBIMOS DATOS
    if (energies_switch.eqv..true.) then
        open(12,file='../results/result_03.dat',status='replace',action='write',iostat=istat)
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
        call velocity_verlet(n_p,x_vector,y_vector,z_vector,&
        vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force_x,force_y,force_z)
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
        else if ((movie_switch.eqv..true.).and.(mod(i,100)==0_sp)) then
            index=index+1;call create_movie(index,x_vector,y_vector,z_vector,n_p)
        end if
    end do
    if (energies_switch.eqv..true.) close(12)

    deallocate(x_vector,y_vector,z_vector)
    deallocate(vx_vector,vy_vector,vz_vector)
    deallocate(force_x,force_y,force_z)

    call cpu_time(time_end)
    write(*,*) 'elapsed time = ',time_end-time_start,'[s]'
end program molecular_dynamic_lennard_jones

subroutine create_movie(index,x_vector,y_vector,z_vector,n_p)
    use module_precision
    implicit none
    integer(sp), intent(in) :: index,n_p
    real(dp),    intent(in) :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
    character(len=24)       :: file_name
    character(len=2)        :: index_str
    integer(sp)             :: i,istat
    50 format(2(A12,x),A12);51 format(2(E12.4,x),E12.4)
    write (index_str,'(I2)') index
    file_name='../results/picture'//trim(index_str)//'.dat'
    open(52,file=file_name,status='replace',action='write',iostat=istat)
    if (istat/=0) write(*,*) 'ERROR! istat(52file) = ',istat
    write(52,50) 'rx_fcc','ry_fcc','rz_fcc'
    do i=1,n_p;write(52,51) x_vector(i),y_vector(i),z_vector(i);end do;close(52)
end subroutine create_movie

subroutine statistics_variables

end subroutine