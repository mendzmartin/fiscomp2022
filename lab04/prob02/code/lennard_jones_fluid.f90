! make clean && make lennard_jones_fluid.o && ./lennard_jones_fluid.o
program lennard_jones_fluid
    use module_precision;use module_mc_lennard_jones
    implicit none
    integer(sp), parameter   :: MC_step=10_sp,MC_step_trans=10_sp ! total and transitory Monte Carlo step
    ! VARIABLES y PARAMETROS GENERALES
    integer(sp), parameter   :: n_p=256_sp                             ! cantidad de partículas
    real(dp),    parameter   :: T_adim=1._dp                           ! temperatura de referencia adimensional
    real(dp),    parameter   :: density=0.8_dp                         ! densidad (particulas/volumen)
    real(dp),    parameter   :: r_cutoff=2.5_dp,mass=1._dp             ! radio de corte de interacciones y masa    
    real(dp),    allocatable :: x_vector(:),y_vector(:),z_vector(:)    ! componentes de las posiciones/particula
    real(dp),    allocatable :: x_vector_noPBC(:),y_vector_noPBC(:),&  ! componentes de la posición sin PBC
                                z_vector_noPBC(:)
    ! real(dp),    allocatable :: force_x(:),force_y(:),force_z(:)       ! componentes de la fuerza/particula
    integer(sp)              :: istat                        ! loop index
    real(dp)                 :: time_end,time_start                    ! tiempos de CPU
    ! VARIABLES PARA COMPUTAR ENERGÍA POTENCIAL Y PRESIÓN OSMÓTICA
    real(dp)                 :: U_adim,P_adim                          ! observables

    ! comienzo a medir tiempo de cpu
    call cpu_time(time_start)

    ! alloco memoria e inicializo arreglos en cero.
    allocate(x_vector(n_p),y_vector(n_p),z_vector(n_p))
    x_vector(:)=0._dp;y_vector(:)=0._dp;z_vector(:)=0._dp
    allocate(x_vector_noPBC(n_p),y_vector_noPBC(n_p),z_vector_noPBC(n_p))
    x_vector_noPBC(:)=0._dp;y_vector_noPBC(:)=0._dp;z_vector_noPBC(:)=0._dp
    ! allocate(force_x(n_p),force_y(n_p),force_z(n_p))
    ! force_x(:)=0._dp;force_y(:)=0._dp;force_z(:)=0._dp

    ! generamos configuración inicial (FCC structure)
    call initial_lattice_configuration(n_p,density,x_vector,y_vector,z_vector,2)
    ! cargamos posiciones iniciales sin PBC
    x_vector_noPBC(:)=x_vector(:);y_vector_noPBC(:)=y_vector(:);z_vector_noPBC(:)=z_vector(:)
    ! computamos energía interna, fuerzas y presión en el tiempo inicial
    U_adim=u_lj_total(n_p,x_vector,y_vector,z_vector,r_cutoff,density)
    P_adim=osmotic_pressure(n_p,density,mass,r_cutoff,x_vector,y_vector,z_vector)
    ! call f_lj_total(x_vector,y_vector,z_vector,r_cutoff,n_p,density,force_x,force_y,force_z)

    ! TRANSITORIO
    ! evolucionamos sistema con método de metropolis
    open(10,file='../results/results_naive.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat
    write(10,'(2(A12,x),A12)') 'MC_step','U_adim','P_adim'
    call lennard_jones_relaxation(n_p,MC_step,MC_step_trans,10,&
        x_vector,y_vector,z_vector,x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
        U_adim,T_adim,P_adim,mass,r_cutoff,density)

    ! terminamos de medir tiempo de cpu
    call cpu_time(time_end)

    write(*,*) 'elapsed time = ',(time_end-time_start),'[s]'

    ! liberamos memoria
    deallocate(x_vector,y_vector,z_vector)
    deallocate(x_vector_noPBC,y_vector_noPBC,z_vector_noPBC)
    ! deallocate(force_x,force_y,force_z)
    
end program lennard_jones_fluid

subroutine lennard_jones_relaxation(n_p,MC_step,MC_step_trans,file_num,&
    x_vector,y_vector,z_vector,x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
    U_adim,T_adim,P_adim,mass,r_cutoff,density)

    use module_precision;use module_mc_lennard_jones
    implicit none
    integer(sp), intent(in)     :: n_p,file_num,MC_step,MC_step_trans
    real(dp),    intent(inout)  :: x_vector(n_p),y_vector(n_p),z_vector(n_p)
    real(dp),    intent(inout)  :: x_vector_noPBC(n_p),y_vector_noPBC(n_p),z_vector_noPBC(n_p)
    real(dp),    intent(in)     :: mass,r_cutoff,density,T_adim
    real(dp),    intent(inout)  :: U_adim,P_adim
    integer(sp)                 :: i,j
    real(dp)                    :: U_med_adim,sigma_U,error_U
    real(dp)                    :: P_med_adim,sigma_P,error_P
    real(dp)                    :: s0,s1_U,s2_U,s1_P,s2_P                 ! variables para hacer estadística (momentos)

    90 format(I12,x,E12.4,x,E12.4)
    s0=0._dp;s1_U=0._dp;s2_U=0._dp;s1_P=0._dp;s2_P=0._dp
    do i=1,MC_step
        write(*,*) i
        call evolution_monte_carlo(n_p,x_vector,y_vector,z_vector,&
            x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
            U_adim,T_adim,r_cutoff,density)
        do j=1,n_p
            if (x_vector(j)==0._dp) write(*,*) 'x=0'
        end do
        U_adim=u_lj_total(n_p,x_vector,y_vector,z_vector,r_cutoff,density)
        P_adim=osmotic_pressure(n_p,density,mass,r_cutoff,x_vector,y_vector,z_vector)
        ! datos para hacer estadística en steady state
        if (i>=MC_step_trans) then
            s0=s0+1._dp
            s2_U=s2_U+U_adim*U_adim;s1_U=s1_U+U_adim ! energía
            s2_P=s2_P+P_adim*P_adim;s1_P=s1_P+P_adim ! presión osmótica
        end if
        write(file_num,90) i,U_adim,P_adim
    end do

    !calculamos media,desviacion estándar, varianza y error
    U_med_adim=s1_U*(1._dp/s0)
    sigma_U=sqrt((s2_U-s0*U_med_adim*U_med_adim)*(1._dp/(s0-1._dp)))
    error_U=sigma_U*(1._dp/sqrt(s0-1._dp))

    P_med_adim=s1_P*(1._dp/s0)
    sigma_P=sqrt((s2_P-s0*P_med_adim*P_med_adim)*(1._dp/(s0-1._dp)))
    error_P=sigma_P*(1._dp/sqrt(s0-1._dp))

    write(*,'(A9,I2,A2,I2)')  'Lattice=',n_p,'x',n_p
    write(*,'(2(A14,E12.4))') 'P_med_adim=',P_med_adim,'error_M',error_P
    write(*,'(2(A14,E12.4))') 'U_med_adim=',U_med_adim,'error_U',error_U

end subroutine lennard_jones_relaxation