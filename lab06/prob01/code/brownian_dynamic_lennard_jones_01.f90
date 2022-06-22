! make clean && make brownian_dynamic_lennard_jones_01.o && ./brownian_dynamic_lennard_jones_01.o
program brownian_dynamic_lennard_jones_01
    use module_precision;use module_bd_lennard_jones
    implicit none
    integer(sp), parameter   :: n_p=256_sp                             ! cantidad de partículasa
    real(dp),    parameter   :: delta_time=0.001_dp                    ! paso temporal
    integer(sp), parameter   :: time_eq=2000_sp,&                      ! pasos de equilibración
                                time_run=1000_sp                       ! pasos de evolucion en el estado estacionario
    real(dp),    parameter   :: T_adim_ref=1.1_dp                      ! temperatura de referencia adimensional
    real(dp),    parameter   :: density=0.8_dp                         ! densidad (particulas/volumen)
    real(dp),    parameter   :: r_cutoff=2.5_dp,mass=1._dp             ! radio de corte de interacciones y masa     
    real(dp),    allocatable :: x_vector(:),y_vector(:),z_vector(:)    ! componentes de las posiciones/particula
    real(dp),    allocatable :: force_x(:),force_y(:),force_z(:)       ! componentes de la fuerza/particula
    integer(sp)              :: i,istat,index                          ! loop index
    real(dp)                 :: U_adim,time,press  ! observables
    real(dp)                 :: time_end,time_start                    ! tiempos de CPU
    real(dp),    parameter   :: dinamic_viscosity=2.87_dp,pi=4._dp*atan(1._dp)
    real(dp),    parameter   :: diffusion_coeff=T_adim_ref*(1._dp/(3._dp*pi*dinamic_viscosity))
    logical                  :: movie_switch,fcc_init_switch,&
                                energies_switch

    movie_switch       =.false. ! escribir pelicula con partículas en la caja
    fcc_init_switch    =.false. ! escribir estructura fcc inicial
    energies_switch    =.false. ! escribir energías en el estado estacionario

    call cpu_time(time_start)
    22 format(2(E12.4,x),E12.4);23 format(2(A12,x),A12)

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

    ! computamos fuerzas en el tiempo inicial
    allocate(force_x(n_p),force_y(n_p),force_z(n_p))
    force_x(:)=0._dp;force_y(:)=0._dp;force_z(:)=0._dp
    call f_lj_total(x_vector,y_vector,z_vector,r_cutoff,n_p,density,force_x,force_y,force_z)

    ! TRANSITORIO
    do i=1,time_eq
        write(*,*) i
        call evolution_bd(n_p,x_vector,y_vector,z_vector,&
            delta_time,mass,r_cutoff,density,force_x,force_y,force_z,&
            dinamic_viscosity,diffusion_coeff)
    end do

    ! ESTACIONARIO
    time=0._dp
    ! ESCRIBIMOS DATOS
    if (energies_switch.eqv..true.) then
        open(12,file='../results/result_03.dat',status='replace',action='write',iostat=istat)
        write(*,*) 'istat(12file) = ',istat;write(12,23) 'time','pot_ergy','press'
    end if
    U_adim=u_lj_total(n_p,x_vector,y_vector,z_vector,r_cutoff,density)
    press=osmotic_pressure(n_p,density,mass,r_cutoff,x_vector,y_vector,z_vector)
    ! ESCRIBIMOS DATOS
    if (energies_switch.eqv..true.) then
        write(12,22) time,U_adim*(1._dp/real(n_p,dp)),press
    end if

    do i=1,time_run
        write(*,*) time_eq+i
        call evolution_bd(n_p,x_vector,y_vector,z_vector,&
            delta_time,mass,r_cutoff,density,force_x,force_y,force_z,&
            dinamic_viscosity,diffusion_coeff)

        U_adim=u_lj_total(n_p,x_vector,y_vector,z_vector,r_cutoff,density)
        press=osmotic_pressure(n_p,density,mass,r_cutoff,x_vector,y_vector,z_vector)
        time=real(i,dp)*delta_time

        ! ESCRIBIMOS DATOS
        if (energies_switch.eqv..true.) then
            write(12,22) time,U_adim*(1._dp/real(n_p,dp)),press
        else if ((movie_switch.eqv..true.).and.(mod(i,100)==0_sp)) then
            index=index+1;call create_movie(index,x_vector,y_vector,z_vector,n_p)
        end if
    end do
    if (energies_switch.eqv..true.) close(12)

    deallocate(x_vector,y_vector,z_vector)
    deallocate(force_x,force_y,force_z)

    call cpu_time(time_end)
    write(*,*) 'elapsed time = ',time_end-time_start,'[s]'

end program brownian_dynamic_lennard_jones_01

! subrutina para crear película de partículas en la caja
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