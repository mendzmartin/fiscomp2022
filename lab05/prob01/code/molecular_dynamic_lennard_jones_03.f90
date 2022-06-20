! Inciso c)
! make clean && make molecular_dynamic_lennard_jones_03.o && ./molecular_dynamic_lennard_jones_03.o
program molecular_dynamic_lennard_jones_03
    use module_precision;use module_md_lennard_jones
    implicit none
    integer(sp), parameter   :: n_p=256_sp                             ! cantidad de partículasa
    real(dp),    parameter   :: delta_time=0.005_dp                    ! paso temporal
    integer(sp), parameter   :: time_eq=1000_sp,time_scal=50_sp        ! pasos de equilibración y de escaleo de veloc.
    real(dp),    parameter   :: T_adim_ref=1.1_dp                      ! temperatura de referencia adimensional
    real(dp),    parameter   :: density=0.8_dp                         ! densidad (particulas/volumen)
    real(dp),    parameter   :: r_cutoff=2.5_dp,mass=1._dp             ! radio de corte de interacciones,masa     
    real(dp),    allocatable :: x_vector(:),y_vector(:),z_vector(:)    ! componentes de las posiciones/particula
    real(dp),    allocatable :: vx_vector(:),vy_vector(:),vz_vector(:) ! componentes de la velocidad/particula
    real(dp),    allocatable :: vtot_vector(:)
    real(dp),    allocatable :: force_x(:),force_y(:),force_z(:)       ! componentes de la fuerza/particula
    integer(sp)              :: i,index,istat                          ! loop index
    real(dp)                 :: T_adim
    real(dp)                 :: vx_mc,vy_mc,vz_mc                      ! componentes de la velocidad del centro de masas
    real(dp)                 :: time_end,time_start                    ! tiempos de CPU

    real(dp),   allocatable  :: probability_vx(:),probability_vy(:),&
                                probability_vz(:),probability_vtot(:)
    real(dp),   allocatable  :: exact_probability_vx(:),exact_probability_vy(:),&
                                exact_probability_vz(:),exact_probability_vtot(:)
    real(dp),   allocatable  :: variable_vx(:),variable_vy(:),&
                                variable_vz(:),variable_vtot(:)
    integer(sp)              :: n_bins   ! numbers of bins

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

    index=0
    ! TRANSITORIO
    do i=1,time_eq
        index=index+1
        write(*,*) 'paso temporal =',index,' de',time_eq
        if (mod(i,time_scal)==0_sp) call rescaling_velocities(n_p,vx_vector,vy_vector,vz_vector,T_adim_ref,mass)
        call velocity_verlet(n_p,x_vector,y_vector,z_vector,&
        vx_vector,vy_vector,vz_vector,delta_time,mass,r_cutoff,density,force_x,force_y,force_z)
        ! velocity center of mass to zero
        vx_mc=sum(vx_vector(:))*(1._dp/real(n_p,dp));vx_vector(:)=(vx_vector(:)-vx_mc)
        vy_mc=sum(vy_vector(:))*(1._dp/real(n_p,dp));vy_vector(:)=(vy_vector(:)-vy_mc)
        vz_mc=sum(vz_vector(:))*(1._dp/real(n_p,dp));vz_vector(:)=(vz_vector(:)-vz_mc)
        T_adim=temperature(n_p,mass,vx_vector,vy_vector,vz_vector)
    end do

    deallocate(x_vector,y_vector,z_vector)
    deallocate(force_x,force_y,force_z)

    ! computamos distribuciones exactas de velocidad (Maxwell-Boltzmann y Gauss)
    allocate(vtot_vector(n_p));vtot_vector(:)=0._dp
    allocate(exact_probability_vx(n_p),exact_probability_vy(n_p),exact_probability_vz(n_p))
    allocate(exact_probability_vtot(n_p))
    open(10,file='../results/exact_velocities_distributions.dat',status='replace',action='write',iostat=istat)
    if (istat/=0) write(*,*) 'ERROR! istat(11file) = ',istat
    24 format(7(E12.4,x),E12.4);25 format(7(A12,x),A12)
    write(10,25) 'vx','p(vx)','vy','p(vy)','vz','p(vz)','v_tot','p(vtot)'
    do i=1,n_p
        exact_probability_vx(i)=sqrt(mass*0.125_dp*(1._dp/atan(1._dp))*(1._dp/T_adim_ref))*&
            exp(-mass*vx_vector(i)*vx_vector(i)*0.5_dp*(1._dp/T_adim_ref))

        exact_probability_vy(i)=sqrt(mass*0.125_dp*(1._dp/atan(1._dp))*(1._dp/T_adim_ref))*&
            exp(-mass*vy_vector(i)*vy_vector(i)*0.5_dp*(1._dp/T_adim_ref))

        exact_probability_vz(i)=sqrt(mass*0.125_dp*(1._dp/atan(1._dp))*(1._dp/T_adim_ref))*&
            exp(-mass*vz_vector(i)*vz_vector(i)*0.5_dp*(1._dp/T_adim_ref))

        vtot_vector(i)=sqrt(vx_vector(i)*vx_vector(i)+vy_vector(i)*vy_vector(i)+vz_vector(i)*vz_vector(i))

        exact_probability_vtot(i)=sqrt(0.5_dp*(1._dp/atan(1._dp))*(mass*(1._dp/T_adim_ref))**3)*&
            vtot_vector(i)*vtot_vector(i)*exp(-mass*vtot_vector(i)*vtot_vector(i)*0.5_dp*(1._dp/T_adim_ref))
        
        write(10,24) vx_vector(i),exact_probability_vx(i),vy_vector(i),exact_probability_vy(i),&
                     vz_vector(i),exact_probability_vz(i),vtot_vector(i),exact_probability_vtot(i)
    end do
    deallocate(exact_probability_vx,exact_probability_vy,exact_probability_vz)
    deallocate(exact_probability_vtot)
    close(10)

    ! hacemos histogramas de las componentes de la velocidad
    n_bins=50
    allocate(probability_vx(n_bins),variable_vx(n_bins+1))
    call histogram(vx_vector,n_p,variable_vx,probability_vx,n_bins)

    allocate(probability_vy(n_bins),variable_vy(n_bins+1))
    call histogram(vy_vector,n_p,variable_vy,probability_vy,n_bins)

    allocate(probability_vz(n_bins),variable_vz(n_bins+1))
    call histogram(vz_vector,n_p,variable_vz,probability_vz,n_bins)

    open(10,file='../results/components_velocities_histogram.dat',status='replace',action='write',iostat=istat)
    if (istat/=0) write(*,*) 'ERROR! istat(11file) = ',istat
    20 format(5(E12.4,x),E12.4);21 format(5(A12,x),A12)
    write(10,21) 'vx','p(vx)','vy','p(vy)','vz','p(vz)'
    do i = 1,n_bins
        write(10,20) variable_vx(i),probability_vx(i)*(1._dp/(variable_vx(i+1)-variable_vx(i))),&
                     variable_vy(i),probability_vy(i)*(1._dp/(variable_vy(i+1)-variable_vy(i))),&
                     variable_vz(i),probability_vz(i)*(1._dp/(variable_vz(i+1)-variable_vz(i)))
        ! write(10,20) variable_vx(i),probability_vx(i),&
        !              variable_vy(i),probability_vy(i),&
        !              variable_vz(i),probability_vz(i)
    end do
    deallocate(vx_vector,vy_vector,vz_vector)
    deallocate(probability_vx,variable_vx)
    deallocate(probability_vy,variable_vy)
    deallocate(probability_vz,variable_vz)

    ! hacemos histograma de la velocidad total
    n_bins=50
    allocate(probability_vtot(n_bins),variable_vtot(n_bins+1))
    call histogram(vtot_vector,n_p,variable_vtot,probability_vtot,n_bins)

    open(11,file='../results/total_velocities_histogram.dat',status='replace',action='write',iostat=istat)
    if (istat/=0) write(*,*) 'ERROR! istat(11file) = ',istat
    22 format(E12.4,x,E12.4);23 format(A12,x,A12)
    write(11,23) 'vtot','p(vtot)'
    do i = 1,n_bins
        write(11,22) variable_vtot(i),probability_vtot(i)*(1._dp/(variable_vtot(i+1)-variable_vtot(i)))
    end do
    close(11)
    deallocate(vtot_vector)
    deallocate(probability_vtot,variable_vtot)

    call cpu_time(time_end)
    write(*,*) 'elapsed time = ',time_end-time_start,'[s]'
end program molecular_dynamic_lennard_jones_03

! subrutina para crear e imprimir histograma
subroutine histogram(x_vector,x_dim,variable,probability,n_bins)
    use module_precision
    implicit none
    integer(sp),      intent(in)    :: x_dim           ! dimension
    integer(sp),      intent(in)    :: n_bins
    real(dp),         intent(in)    :: x_vector(x_dim) ! data to do histogram
    real(dp),         intent(inout) :: variable(n_bins+1),probability(n_bins)

    integer(sp), allocatable :: counter(:)      ! counter vector of bins
    real(sp)                 :: max_value       ! maximun counter value
    real(dp)                 :: min_bin_point,max_bin_point
    real(dp)                 :: bins_step       ! step of points between bins
    integer(sp)              :: i,j       ! loop and control variables

    allocate(counter(n_bins))

    ! armamos el vector de bins (o variable(:)) en el rango [x_min,x_max]
    max_bin_point=maxval(x_vector(:));min_bin_point=minval(x_vector(:))
    bins_step=abs(max_bin_point-min_bin_point)*(1._dp/n_bins)
    do i=1,n_bins+1;variable(i)=min_bin_point+bins_step*real(i-1,dp);end do

    ! llenamos el vector contador de bins
    do i=1,n_bins
        counter(i)=0
        do j=1,x_dim
            if ((variable(i)<=x_vector(j)).and.(variable(i+1)>=x_vector(j))) then
                counter(i)=counter(i)+1
            endif
    enddo;enddo

    max_value=abs(real(maxval(counter(:)),dp))
    if (max_value==0._dp) write(*,*) 'math error'

    ! escribimos distribución de probabilidad
    probability(:)=real(counter(:),dp)*(1._dp/max_value)

    deallocate(counter)
end subroutine histogram