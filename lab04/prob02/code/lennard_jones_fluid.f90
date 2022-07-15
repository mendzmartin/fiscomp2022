! make clean && make lennard_jones_fluid.o && ./lennard_jones_fluid.o
program lennard_jones_fluid
    use module_precision;use module_mc_lennard_jones
    implicit none
    integer(sp), parameter   :: MC_step=10000_sp,MC_step_trans=3000_sp ! Monte Carlo step total and transitory
    ! VARIABLES y PARAMETROS GENERALES
    integer(sp), parameter   :: n_p=256_sp                             ! cantidad de partículasa
    real(dp),    parameter   :: delta_time=0.001_dp                    ! paso temporal
    integer(sp), parameter   :: time_eq=3000_sp,&                     ! pasos de equilibración
                                time_run=10000_sp,&                    ! pasos de evolucion en el estado estacionario
                                ensamble_step=10_sp                    ! pasos de evolución para promedio en ensamble
    real(dp),    parameter   :: T_adim_ref=1._dp                       ! temperatura de referencia adimensional
    real(dp),    parameter   :: density=0.8_dp                         ! densidad (particulas/volumen)
    real(dp),    parameter   :: r_cutoff=2.5_dp,mass=1._dp             ! radio de corte de interacciones y masa    
    real(dp),    parameter   :: dinamic_viscosity=2.87_dp
    real(dp),    parameter   :: pi=4._dp*atan(1._dp)
    real(dp),    parameter   :: diffusion_coeff=T_adim_ref*&
                                 (1._dp/(3._dp*pi*dinamic_viscosity)) 
    real(dp),    allocatable :: x_vector(:),y_vector(:),z_vector(:)    ! componentes de las posiciones/particula
    real(dp),    allocatable :: x_vector_noPBC(:),y_vector_noPBC(:),&  ! componentes de la posición sin PBC
                                z_vector_noPBC(:)
    real(dp),    allocatable :: force_x(:),force_y(:),force_z(:)       ! componentes de la fuerza/particula
    integer(sp)              :: i,j,istat,index                        ! loop index
    real(dp)                 :: time_end,time_start                    ! tiempos de CPU
    ! VARIABLES PARA COMPUTAR ENERGÍA POTENCIAL
    real(dp)                 :: U_adim,time,press  ! observables
    real(dp)                 :: U_med,var_U,err_U
    real(dp)                 :: s1_U,s2_U
    ! VARIABLES PARA COMPUTAR PRESIÓN OSMÓTICA
    real(dp)                 :: press_med,var_press,err_press
    real(dp)                 :: s1_press,s2_press

    call cpu_time(time_start)

    allocate(x_vector(n_p),y_vector(n_p),z_vector(n_p))
    x_vector(:)=0._dp;y_vector(:)=0._dp;z_vector(:)=0._dp
    allocate(x_vector_noPBC(n_p),y_vector_noPBC(n_p),z_vector_noPBC(n_p))
    x_vector_noPBC(:)=0._dp;y_vector_noPBC(:)=0._dp;z_vector_noPBC(:)=0._dp
    allocate(force_x(n_p),force_y(n_p),force_z(n_p))
    force_x(:)=0._dp;force_y(:)=0._dp;force_z(:)=0._dp

    ! generamos configuración inicial (FCC structure)
    call initial_lattice_configuration(n_p,density,x_vector,y_vector,z_vector,2)
    x_vector_noPBC(:)=x_vector(:);y_vector_noPBC(:)=y_vector(:);z_vector_noPBC(:)=z_vector(:)
    ! computamos fuerzas en el tiempo inicial
    call f_lj_total(x_vector,y_vector,z_vector,r_cutoff,n_p,density,force_x,force_y,force_z)

    ! TRANSITORIO
    do i=1,time_eq
        call evolution_monte_carlo(n_p,x_vector,y_vector,z_vector,&
        x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
        U_adim,T_adim,r_cutoff,density)
        call evolution_bd(n_p,x_vector,y_vector,z_vector,&
            x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
            delta_time,mass,r_cutoff,density,force_x,force_y,force_z,&
            dinamic_viscosity,diffusion_coeff)
    end do

    
end program lennard_jones_fluid

subroutine lennard_jones_relaxation(n,MC_step,MC_step_trans,m,file_num,aux_matrix_pbc,T_start,T_end,Tc_adim,U_adim)
    use module_precision;use module_mc_lennard_jones
    implicit none 
    real(dp),    intent(in)     :: T_end,T_start,Tc_adim
    integer(sp), intent(in)     :: n,m,file_num,MC_step,MC_step_trans
    integer(sp), intent(inout)  :: aux_matrix_pbc(n+2,n+2)
    real(dp),    intent(inout)  :: U_adim,Madim
    integer(sp)                 :: i,j
    real(dp)                    :: T_step,T_adim
    real(dp)                    :: U_med_adim,sigma_U,error_U
    real(dp)                    :: M_med_adim,sigma_M,error_M
    real(dp)                    :: s0,s1_U,s2_U,s1_M,s2_M ! variables para hacer estadística

    20 format(I12,x,E12.4,x,E12.4)
    T_step=abs(T_end-T_start)*(1._dp/real(m-1_sp,dp))
    ! termalizo hasta una temperatura anterior a la buscada
    do j=1,m-1
        T_adim=T_start+T_step*real(j-1_sp,dp)
        write(*,'(A12,E12.4)') 'T_adim=',T_adim
        do i=1,MC_step
            call evolution_monte_carlo(n_p,x_vector,y_vector,z_vector,&
                x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
                U_adim,T_adim,r_cutoff,density)
        end do
    end do
    ! termalizo a la temperatura buscada
    write(*,'(A12,E12.4)') 'T_adim=',T_end
    s0=0._dp;s1_U=0._dp;s2_U=0._dp;s1_M=0._dp;s2_M=0._dp
    do i=1,MC_step
        call evolution_monte_carlo(n_p,x_vector,y_vector,z_vector,&
            x_vector_noPBC,y_vector_noPBC,z_vector_noPBC,&
            U_adim,T_adim,r_cutoff,density)
        ! datos para hacer estadística en steady state
        if (i>=MC_step_trans) then
            s0=s0+1._dp
            s2_U=s2_U+U_adim*U_adim;s1_U=s1_U+U_adim ! energía
        end if
        write(file_num,20) i,U_adim
    end do

    !calculamos media,desviacion estándar, varianza y error
    U_med_adim=s1_U*(1._dp/s0)
    sigma_U=sqrt((s2_U-s0*U_med_adim*U_med_adim)*(1._dp/(s0-1._dp)))
    error_U=sigma_U*(1._dp/sqrt(s0-1._dp))
    M_med_adim=s1_M*(1._dp/s0)
    sigma_M=sqrt((s2_M-s0*M_med_adim*M_med_adim)*(1._dp/(s0-1._dp)))
    error_M=sigma_M*(1._dp/sqrt(s0-1._dp))

    write(*,'(A9,I2,A2,I2)')  'Lattice=',n,'x',n
    write(*,'(A14,E12.4)')    'M_exact_adim=',M_exact_adim(n,T_end,Tc_adim)
    write(*,'(2(A14,E12.4))') 'M_med_adim=',M_med_adim,'error_M',error_M
    write(*,'(2(A14,E12.4))') 'U_med_adim=',U_med_adim,'error_U',error_U

end subroutine lennard_jones_relaxation