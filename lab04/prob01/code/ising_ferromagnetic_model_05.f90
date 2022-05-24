! P01.d (calculamos los cumulantes de binder)
! make clean && make ising_ferromagnetic_model_05.o && ./ising_ferromagnetic_model_05.o
program ising_ferromagnetic_model_05
    use module_precision;use module_2D_ferromagnetic_ising
    implicit none

    integer(sp), parameter   :: n=10_sp                                 ! sitios de red por dimension
    integer(sp), parameter   :: MC_step=5000_sp,MC_step_trans=3000_sp  ! Monte Carlo step total and transitory
    integer(sp), parameter   :: m1=25_sp,m2=50_sp                    ! puntos p/ mallado fino y grueso
    real(dp),    parameter   :: Tmin_adim=0._dp,Tmax_adim=10._dp        ! temperatura adimensional
    real(dp),    parameter   :: Tc_adim=2.2692_dp                       ! temperatura de Curie adimensional
    real(dp),    parameter   :: deltaT_adim=Tc_adim*0.5_dp            ! intervalo para incremetar ptos
    integer(sp), allocatable :: aux_matrix_pbc(:,:)
    integer(sp)              :: i,istat
    real(dp)                 :: U_adim,U_med_adim,sigma_U,error_U       ! Energía interna
    real(dp)                 :: Madim,M_med_adim,sigma_M,error_M,Mexact ! Magenitación
    real(dp)                 :: susc_adim                               ! Susceptibilidad magnética
    real(dp)                 :: cv                                      ! calor específico 
    real(dp)                 :: T_adim,T_step1,T_step2
    real(dp)                 :: time_start,time_end                     ! tiempos de CPU
    real(dp)                 :: binder_cumulant

    open(10,file='../results/result_01d_binder_cumulant_10x10.dat',status='replace',action='write',iostat=istat)
    !open(10,file='../results/result_01d_binder_cumulant_20x20.dat',status='replace',action='write',iostat=istat)
    !open(10,file='../results/result_01d_binder_cumulant_40x40.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat
    21 format(A12,x,A12);write(10,21) 'T/Tc','binder_cumulant'
    20 format(E12.4,x,E12.4)

    call cpu_time(time_start)

    allocate(aux_matrix_pbc(n+2,n+2))
    ! genero configuracion inicial (random,descorrelacionada)
    call initial_spins_configuration(1_sp,n,aux_matrix_pbc)
    ! calculamos energía interna (configuración inicial)
    call average_energy(aux_matrix_pbc,n,U_adim)
    ! calculamos magnetización inicial
    Madim=M_adim(aux_matrix_pbc,n)

    ! mallado grueso
    do i=1,m1
        T_step1=abs(Tc_adim-deltaT_adim-Tmin_adim)*(1._dp/real(m1-1,dp))
        T_adim=Tmin_adim+T_step1*real(i-1,dp)
        call rlx_ising(n,aux_matrix_pbc,MC_step,MC_step_trans,T_adim,Tc_adim,&
                       U_adim,U_med_adim,sigma_U,error_U,&
                       Madim,M_med_adim,sigma_M,error_M,Mexact,&
                       susc_adim,cv,binder_cumulant)
        write(10,20) T_adim*(1._dp/Tc_adim),binder_cumulant
    end do
    write(*,*) 'primer bucle terminado...'

    ! mallado fino
    do i=1,m2
        T_step2=2_dp*deltaT_adim*(1._dp/real(m2-1,dp))
        T_adim=(Tc_adim-deltaT_adim)+T_step2*real(i-1,dp)
        call rlx_ising(n,aux_matrix_pbc,MC_step,MC_step_trans,T_adim,Tc_adim,&
                       U_adim,U_med_adim,sigma_U,error_U,&
                       Madim,M_med_adim,sigma_M,error_M,Mexact,&
                       susc_adim,cv,binder_cumulant)
        write(10,20) T_adim*(1._dp/Tc_adim),binder_cumulant
    end do
    write(*,*) 'segundo bucle terminado...'

    ! mallado grueso
    do i=1,m1
        T_step1=abs(Tmax_adim-(Tc_adim+deltaT_adim))*(1._dp/real(m1-1,dp))
        T_adim=(Tc_adim+deltaT_adim)+T_step1*real(i-1,dp)
        call rlx_ising(n,aux_matrix_pbc,MC_step,MC_step_trans,T_adim,Tc_adim,&
                       U_adim,U_med_adim,sigma_U,error_U,&
                       Madim,M_med_adim,sigma_M,error_M,Mexact,&
                       susc_adim,cv,binder_cumulant)
        write(10,20) T_adim*(1._dp/Tc_adim),binder_cumulant
    end do

    close(10)
    deallocate(aux_matrix_pbc)

    call cpu_time(time_end); write(*,*) 'elapsed time = ',(time_end-time_start)
end program ising_ferromagnetic_model_05

subroutine rlx_ising(n,aux_matrix_pbc,MC_step,MC_step_trans,T_adim,Tc_adim,&
                     U_adim,U_med_adim_mom,sigma_U_mom,error_U_mom,&
                     Madim,M_med_adim_mom,sigma_M_mom,error_M_mom,Mexact,&
                     susc_adim,cv,binder_cumulant)
    use module_precision;use module_2D_ferromagnetic_ising
    implicit none

    integer(sp), intent(in)     :: n,MC_step,MC_step_trans
    real(dp),    intent(in)     :: T_adim,Tc_adim
    integer(sp), intent(inout)  :: aux_matrix_pbc(n+2,n+2)
    real(dp),    intent(inout)  :: U_adim,U_med_adim_mom,sigma_U_mom,error_U_mom       ! Energía interna
    real(dp),    intent(inout)  :: Madim,M_med_adim_mom,sigma_M_mom,error_M_mom,Mexact        ! Magenitación
    real(dp),    intent(inout)  :: susc_adim                               ! Susceptibilidad magnética
    real(dp),    intent(inout)  :: cv                                      ! calor específico 
    real(dp),    intent(inout)  :: binder_cumulant

    integer(sp), parameter   :: m_exp=50_sp                ! numero de experimentos
    real(dp)                 :: s0,s1_U,s2_U,s1_M,s2_M,s4_M   ! variables para hacer estadística
    real(dp)                 :: s4_M_mom,s2_M_mom
    integer(sp)              :: i,j
    real(dp),    allocatable :: U_med_adim_vector(:),M_med_adim_vector(:)

    allocate(U_med_adim_vector(m_exp),M_med_adim_vector(m_exp))
    s0=real(MC_step-MC_step_trans,dp)

    s4_M_mom=0._dp;s2_M_mom=0._dp
    do j=1,m_exp
        s1_U=0._dp;s2_U=0._dp;s1_M=0._dp;s2_M=0._dp;s4_M=0._dp
        do i=1,MC_step
            call MC_step_relaxation(1_sp,n,aux_matrix_pbc,T_adim,U_adim)
            Madim=M_adim(aux_matrix_pbc,n) ! Magnetización
            ! datos para hacer estadística en steady state
            if (i>=MC_step_trans) then
                s2_U=s2_U+U_adim*U_adim;s1_U=s1_U+U_adim ! primer y segundo momento energía
                s2_M=s2_M+Madim*Madim;s1_M=s1_M+Madim    ! primer y segundo momento magnetización
                s4_M=s4_M+Madim*Madim*Madim*Madim        ! cuarto momento magnetización
            end if
        end do
        ! calculamos valores medios
        U_med_adim_vector(j)=s1_U*(1._dp/s0)
        M_med_adim_vector(j)=s1_M*(1._dp/s0)
        s4_M_mom=s4_M_mom+s4_M*(1._dp/s0);s2_M_mom=s2_M_mom+s2_M*(1._dp/s0)
    end do

    ! calculamos media de medias,desviacion estándar, varianza y error
    U_med_adim_mom=(1._dp/real(m_exp))*sum(U_med_adim_vector(:))
    U_med_adim_vector(:)=(U_med_adim_vector(:)-U_med_adim_mom)
    U_med_adim_vector(:)=U_med_adim_vector(:)*U_med_adim_vector(:)
    sigma_U_mom=sqrt((1._dp/real(m_exp))*sum(U_med_adim_vector(:)))
    error_U_mom=sigma_U_mom*(1._dp/sqrt(real(m_exp)))

    M_med_adim_mom=(1._dp/real(m_exp))*sum(M_med_adim_vector(:))
    M_med_adim_vector(:)=(M_med_adim_vector(:)-M_med_adim_mom)
    M_med_adim_vector(:)=M_med_adim_vector(:)*M_med_adim_vector(:)
    sigma_M_mom=sqrt((1._dp/real(m_exp))*sum(M_med_adim_vector(:)))
    error_M_mom=sigma_M_mom*(1._dp/sqrt(real(m_exp)))
    if (T_adim==0._dp) then;cv=0._dp;susc_adim=0._dp
    else
        cv=(1._dp/T_adim)*sigma_U_mom*sigma_U_mom
        susc_adim=(1._dp/T_adim)*sigma_M_mom*sigma_M_mom
    end if
    Mexact=M_exact_adim(n,T_adim,Tc_adim)
    binder_cumulant=1._dp-real(m_exp,dp)*s4_M_mom*(1._dp/(3._dp*s2_M_mom*s2_M_mom))

    deallocate(U_med_adim_vector,M_med_adim_vector)

end subroutine rlx_ising