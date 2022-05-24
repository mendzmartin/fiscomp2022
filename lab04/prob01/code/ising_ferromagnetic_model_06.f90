! P01.e (calculamos función de autocorrelación de la energía y magnetización)
! make clean && make ising_ferromagnetic_model_06.o && ./ising_ferromagnetic_model_06.o
program ising_ferromagnetic_model_06
    use module_precision;use module_2D_ferromagnetic_ising
    implicit none

    integer(sp), parameter   :: n=10_sp                                 ! sitios de red por dimension
    integer(sp), parameter   :: MC_step_trans=10000_sp                   ! Monte Carlo step transitory
    integer(sp), parameter   :: m=100_sp                                ! puntos p/ mallado fino
    integer(sp), parameter   :: tau_corr=10_sp                        ! tiempo máximo de correlación
    integer(sp), parameter   :: num_of_terms=1000_sp                   ! cantidad de terminos para la última autocorr
    integer(sp), parameter   :: MC_step=num_of_terms+tau_corr+MC_step_trans ! Monte Carlo step total
    real(dp),    parameter   :: Tmin_adim=0._dp,Tmax_adim=2.2_dp        ! temperatura adimensional
    real(dp),    parameter   :: Tc_adim=2.2692_dp                       ! temperatura de Curie adimensional
    integer(sp), allocatable :: aux_matrix_pbc(:,:)
    real(dp),    allocatable :: autocor_vector_mom(:)                       ! autocorrelación
    integer(sp)              :: i,istat
    real(dp)                 :: U_adim,U_med_adim,sigma_U,error_U       ! Energía interna
    real(dp)                 :: Madim,M_med_adim,sigma_M,error_M,Mexact ! Magenitación
    real(dp)                 :: susc_adim                               ! Susceptibilidad magnética
    real(dp)                 :: cv                                      ! calor específico
    real(dp)                 :: binder_cumulant                         ! cumulante de Binder
    real(dp)                 :: T_adim,T_step
    real(dp)                 :: time_start,time_end                     ! tiempos de CPU

    open(10,file='../results/result_01e_10x10_autocorr_T2.0.dat',status='replace',action='write',iostat=istat)
    ! open(10,file='../results/result_01e_10x10_autocorr_T2.22.dat',status='replace',action='write',iostat=istat)
    ! open(10,file='../results/result_01e_10x10_autocorr_T2.2676.dat',status='replace',action='write',iostat=istat)
    ! open(10,file='../results/result_01e_10x10_autocorr_T2.5.dat',status='replace',action='write',iostat=istat)
    ! open(10,file='../results/result_01e_10x10_autocorr_T3.3.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat
    21 format(A12,x,A12);write(10,21) 'tau_corr','autocorr'
    20 format(I12,x,E12.4)

    call cpu_time(time_start)

    allocate(aux_matrix_pbc(n+2,n+2));allocate(autocor_vector_mom(tau_corr))
    ! genero configuracion inicial (random,descorrelacionada)
    call initial_spins_configuration(1_sp,n,aux_matrix_pbc)
    ! calculamos energía interna (configuración inicial)
    call average_energy(aux_matrix_pbc,n,U_adim)
    ! calculamos magnetización inicial
    Madim=M_adim(aux_matrix_pbc,n)

    do i=1,m
        T_step=abs(Tmax_adim-Tmin_adim)*(1._dp/real(m-1,dp))
        T_adim=Tmin_adim+T_step*real(i-1,dp)
        call rlx_ising(n,aux_matrix_pbc,MC_step,MC_step_trans,T_adim,Tc_adim,Tmax_adim,&
                       U_adim,U_med_adim,sigma_U,error_U,&
                       Madim,M_med_adim,sigma_M,error_M,Mexact,&
                       susc_adim,cv,binder_cumulant,autocor_vector_mom,tau_corr,num_of_terms)
    end do

    do i=1,tau_corr
        write(10,20) i,autocor_vector_mom(i)
    end do

    close(10)
    deallocate(aux_matrix_pbc);deallocate(autocor_vector_mom)

    call cpu_time(time_end); write(*,*) 'elapsed time = ',(time_end-time_start)
end program ising_ferromagnetic_model_06

subroutine rlx_ising(n,aux_matrix_pbc,MC_step,MC_step_trans,T_adim,Tc_adim,Tmax_adim,&
                     U_adim,U_med_adim_mom,sigma_U_mom,error_U_mom,&
                     Madim,M_med_adim_mom,sigma_M_mom,error_M_mom,Mexact,&
                     susc_adim,cv,binder_cumulant,autocor_vector_mom,tau_corr,num_of_terms)
    use module_precision;use module_2D_ferromagnetic_ising
    implicit none

    integer(sp), intent(in)     :: n,MC_step,MC_step_trans,tau_corr,num_of_terms
    real(dp),    intent(in)     :: T_adim,Tc_adim,Tmax_adim
    integer(sp), intent(inout)  :: aux_matrix_pbc(n+2,n+2)
    real(dp),    intent(inout)  :: U_adim,U_med_adim_mom,sigma_U_mom,error_U_mom       ! Energía interna
    real(dp),    intent(inout)  :: Madim,M_med_adim_mom,sigma_M_mom,error_M_mom,Mexact ! Magenitación
    real(dp),    intent(inout)  :: susc_adim                                           ! Susceptibilidad magnética
    real(dp),    intent(inout)  :: cv                                                  ! calor específico 
    real(dp),    intent(inout)  :: binder_cumulant
    real(dp),    intent(inout)  :: autocor_vector_mom(tau_corr)

    integer(sp), parameter   :: m_exp=10_sp                   ! numero de experimentos
    real(dp)                 :: s0,s1_U,s2_U,s1_M,s2_M,s4_M   ! variables para hacer estadística
    real(dp)                 :: s4_M_mom,s2_M_mom,s2_U_mom,sigma_U_mom_v2
    integer(sp)              :: i,j
    real(dp),    allocatable :: U_med_adim_vector(:),M_med_adim_vector(:)
    real(dp),    allocatable :: autocor_vector(:),aux_vector1(:),aux_vector2(:),mask_vector(:)

    allocate(U_med_adim_vector(m_exp),M_med_adim_vector(m_exp))
    allocate(autocor_vector(tau_corr),aux_vector1(tau_corr),aux_vector2(tau_corr),mask_vector(tau_corr))
    s0=real(MC_step-MC_step_trans,dp)

    if (T_adim==Tmax_adim) autocor_vector_mom(:)=0._dp

    s4_M_mom=0._dp;s2_M_mom=0._dp;s2_U_mom=0._dp;sigma_U_mom_v2=0._dp
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
                if (T_adim==Tmax_adim) then
                    call func_autocor(U_adim,i,tau_corr,num_of_terms,&
                                      autocor_vector,aux_vector1,mask_vector,aux_vector2)
                end if
            end if
        end do
        ! calculamos valores medios
        U_med_adim_vector(j)=s1_U*(1._dp/s0)
        M_med_adim_vector(j)=s1_M*(1._dp/s0)
        s4_M_mom=s4_M_mom+s4_M*(1._dp/s0);s2_M_mom=s2_M_mom+s2_M*(1._dp/s0)
        s2_U_mom=s2_U_mom+s2_U*(1._dp/s0)

        !autocorrelación
        autocor_vector(:)=(autocor_vector(:)-U_med_adim_vector(j)*U_med_adim_vector(j))*&
                           (1._dp/(s2_U*(1._dp/s0)-U_med_adim_vector(j)*U_med_adim_vector(j)))
        autocor_vector_mom(:)=autocor_vector_mom(:)+autocor_vector(:)
    end do

    ! calculamos media de medias,desviacion estándar, varianza y error
    U_med_adim_mom=(1._dp/real(m_exp,dp))*sum(U_med_adim_vector(:))
    sigma_U_mom_v2=sqrt(real(m_exp,dp)*s2_U_mom-sum(U_med_adim_vector(:))*sum(U_med_adim_vector(:))*&
                   (1._dp/(real(m_exp,dp)*real(m_exp-1,dp))))
    U_med_adim_vector(:)=(U_med_adim_vector(:)-U_med_adim_mom)
    U_med_adim_vector(:)=U_med_adim_vector(:)*U_med_adim_vector(:)
    sigma_U_mom=sqrt((1._dp/real(m_exp,dp))*sum(U_med_adim_vector(:)))
    error_U_mom=sigma_U_mom*(1._dp/sqrt(real(m_exp,dp)))

    autocor_vector_mom(:)=(1._dp/real(m_exp,dp))*autocor_vector_mom(:)

    M_med_adim_mom=(1._dp/real(m_exp,dp))*sum(M_med_adim_vector(:))
    M_med_adim_vector(:)=(M_med_adim_vector(:)-M_med_adim_mom)
    M_med_adim_vector(:)=M_med_adim_vector(:)*M_med_adim_vector(:)
    sigma_M_mom=sqrt((1._dp/real(m_exp,dp))*sum(M_med_adim_vector(:)))
    error_M_mom=sigma_M_mom*(1._dp/sqrt(real(m_exp,dp)))
    if (T_adim==0._dp) then;cv=0._dp;susc_adim=0._dp
    else
        ! cv=(1._dp/T_adim)*sigma_U_mom*sigma_U_mom
        cv=(1._dp/T_adim)*sigma_U_mom_v2*sigma_U_mom_v2
        susc_adim=(1._dp/T_adim)*sigma_M_mom*sigma_M_mom
    end if
    Mexact=M_exact_adim(n,T_adim,Tc_adim)
    binder_cumulant=1._dp-real(m_exp,dp)*s4_M_mom*(1._dp/(3._dp*s2_M_mom*s2_M_mom))

    deallocate(U_med_adim_vector,M_med_adim_vector)
    deallocate(autocor_vector,aux_vector1,aux_vector2,mask_vector)

end subroutine rlx_ising