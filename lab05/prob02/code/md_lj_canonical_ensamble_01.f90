! make clean && make md_lj_canonical_ensamble_01.o && ./md_lj_canonical_ensamble_01.o
program md_lj_canonical_ensamble_01
    use module_precision;use module_md_lennard_jones
    implicit none
    integer(sp), parameter   :: n_p=500_sp                             ! cantidad de partículasa
    real(dp),    parameter   :: delta_time=0.005_dp                    ! paso temporal
    integer(sp), parameter   :: time_eq=1000_sp,&                      ! pasos de equilibración
                                time_run=1000_sp                       ! pasos de evolucion en el estado estacionario
    real(dp),    parameter   :: T_adim_ref=1.0_dp                      ! temperatura de referencia adimensional
    real(dp),    parameter   :: density=0.8_dp                         ! densidad (particulas/volumen)
    real(dp),    parameter   :: r_cutoff=2.5_dp,mass=1._dp             ! radio de corte de interacciones y masa     
    real(dp),    allocatable :: x_vector(:),y_vector(:),z_vector(:)    ! componentes de las posiciones/particula
    real(dp),    allocatable :: vx_vector(:),vy_vector(:),vz_vector(:) ! componentes de la velocidad/particula
    real(dp),    allocatable :: force_x(:),force_y(:),force_z(:)       ! componentes de la fuerza/particula
    integer(sp)              :: i,istat,index                          ! loop index
    real(dp)                 :: T_adim  ! observables
    real(dp)                 :: vx_mc,vy_mc,vz_mc                      ! componentes de la velocidad del centro de masas
    real(dp)                 :: time,time_end,time_start                    ! tiempos de CPU

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

    ! TRANSITORIO
    do i=1,time_eq
        write(*,*) i
        call rescaling_velocities(n_p,vx_vector,vy_vector,vz_vector,T_adim_ref,mass)
        call velocity_verlet(n_p,x_vector,y_vector,z_vector,&
        vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force_x,force_y,force_z)
        ! velocity center of mass to zero
        vx_mc=sum(vx_vector(:))*(1._dp/real(n_p,dp));vx_vector(:)=(vx_vector(:)-vx_mc)
        vy_mc=sum(vy_vector(:))*(1._dp/real(n_p,dp));vy_vector(:)=(vy_vector(:)-vy_mc)
        vz_mc=sum(vz_vector(:))*(1._dp/real(n_p,dp));vz_vector(:)=(vz_vector(:)-vz_mc)
        T_adim=temperature(n_p,mass,vx_vector,vy_vector,vz_vector)
    end do

    ! ESTACIONARIO
    time=0._dp

    do i=1,time_run
        write(*,*) time_eq+i
        call rescaling_velocities(n_p,vx_vector,vy_vector,vz_vector,T_adim_ref,mass)
        call velocity_verlet(n_p,x_vector,y_vector,z_vector,&
        vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force_x,force_y,force_z)
        time=real(i,dp)*delta_time
    end do

    deallocate(x_vector,y_vector,z_vector)
    deallocate(vx_vector,vy_vector,vz_vector)
    deallocate(force_x,force_y,force_z)

    call cpu_time(time_end)
    write(*,*) 'elapsed time = ',time_end-time_start,'[s]'
end program md_lj_canonical_ensamble_01

subroutine radial_ditribution_function(type_switch,n_p,density,&
    x_vector,y_vector,z_vector,n_bins,g)
    use module_precision;use module_md_lennard_jones
    implicit none
    integer(sp), intent(in) :: type_switch,n_p,n_bins
    real(dp),    intent(in) :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
    real(dp),    intent(in) :: density
    real(dp),    intent(inout) :: g(n_bins) ! funcion distribución

    integer(sp) :: ngr ! contador de particulas
    real(dp)    :: L,rij_pow02,rij,step_bins
    L=(real(n_p,dp)*(1._dp/density))**(1._dp/3._dp)

    if (type_switch==0) then
        ngr=0_sp
        step_bins=L/(2*n_bins)
        g(:)=0._dp
    else if (type_switch==1) then
        ngr=ngr+1_sp
        do j=2,n_p
            do i=1,j-1
                rij_pow02=rel_pos_correction(x_vector(i),y_vector(i),z_vector(i),&
                x_vector(j),y_vector(j),z_vector(j),n_p,density)
                rij=sqrt(rij_pow02)
                if (rij<=L/2) then
                    ig=int(r/step_bins)
                    g(ig)=g(ig)+2._dp
                end if
            end do
        end do
    else if (type_switch==2) then
        do i=1,n_bins
            rij=step_bins*(i+0.5_dp)
            vb=((i+1)**3-i**3)*step_bins**3
            nid=(4/3)*4._dp*atan(1._dp)*vb*density
            g(i)=g(i)/real(ngr*n_p*nid,dp)
        end do
    end if
end subroutine radial_ditribution_function