! P01.e (calculamos función de autocorrelación de la energía y magnetización)
! make clean && make ising_ferromagnetic_model_06.o && ./ising_ferromagnetic_model_06.o
program ising_ferromagnetic_model_06
    use module_precision;use module_2D_ferromagnetic_ising
    implicit none

    integer(sp), parameter   :: n=10_sp                 ! sitios de red por dimension
    integer(sp), parameter   :: MC_step_trans=10000_sp  ! Monte Carlo step transitory
    integer(sp), parameter   :: m=30_sp                 ! puntos p/ para deltas de temperaturas
    integer(sp), parameter   :: tau_corr=1000_sp        ! tiempo máximo de correlación
    integer(sp), parameter   :: num_of_terms=1000000_sp ! cantidad de terminos para la última autocorr
    integer(sp), parameter   :: MC_step=num_of_terms+tau_corr+MC_step_trans-1 ! total Monte Carlo step
    real(dp),    parameter   :: Tmin_adim=0._dp,Tmax_adim=2.5_dp ! temperatura adimensional
    integer(sp), allocatable :: aux_matrix_pbc(:,:)
    real(dp),    allocatable :: autocor_vector_mom_U(:),autocor_vector_mom_M(:)   ! autocorrelación
    integer(sp)              :: i,istat
    real(dp)                 :: U_adim,Madim    ! Energía interna y Magenitación
    real(dp)                 :: T_adim,T_step

    !open(10,file='../results/result_01e_20x20_autocorr_T2.0.dat',status='replace',action='write',iostat=istat)
    !open(10,file='../results/result_01e_10x10_autocorr_T2.22.dat',status='replace',action='write',iostat=istat)
    !open(10,file='../results/result_01e_40x40_autocorr_T2.2676.dat',status='replace',action='write',iostat=istat)
    open(10,file='../results/result_01e_10x10_autocorr_T2.5.dat',status='replace',action='write',iostat=istat)
    !open(10,file='../results/result_01e_20x20_autocorr_T3.3.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat
    21 format(2(A12,x),A12);write(10,21) 'tau_corr','autocorr_U','autocorr_M'
    20 format(I12,x,E12.4,x,E12.4)

    allocate(aux_matrix_pbc(n+2,n+2))
    allocate(autocor_vector_mom_U(tau_corr),autocor_vector_mom_M(tau_corr))
    ! genero configuracion inicial (random,descorrelacionada)
    call initial_spins_configuration(1_sp,n,aux_matrix_pbc)
    ! calculamos energía interna (configuración inicial)
    call average_energy(aux_matrix_pbc,n,U_adim)
    ! calculamos magnetización inicial
    Madim=M_adim(aux_matrix_pbc,n)

    do i=1,m
        T_step=abs(Tmax_adim-Tmin_adim)*(1._dp/real(m-1,dp))
        T_adim=Tmin_adim+T_step*real(i-1,dp)
        if (i==m) T_adim=Tmax_adim
        call rlx_ising(n,aux_matrix_pbc,MC_step,MC_step_trans,T_adim,Tmax_adim,&
                       U_adim,Madim,autocor_vector_mom_U,autocor_vector_mom_M,&
                       tau_corr,num_of_terms)
        write(*,'(I3,x,A3,x,I3)') i,'de',m
    end do

    do i=1,tau_corr
        write(10,20) i,autocor_vector_mom_U(i),autocor_vector_mom_M(i)
    end do

    close(10)
    deallocate(aux_matrix_pbc)
    deallocate(autocor_vector_mom_U);deallocate(autocor_vector_mom_M)
end program ising_ferromagnetic_model_06

subroutine rlx_ising(n,aux_matrix_pbc,MC_step,MC_step_trans,T_adim,Tmax_adim,&
                     U_adim,Madim,autocor_vector_mom_U,autocor_vector_mom_M,&
                     tau_corr,num_of_terms)
    use module_precision;use module_2D_ferromagnetic_ising
    implicit none

    integer(sp), intent(in)     :: n,MC_step,MC_step_trans,tau_corr,num_of_terms
    real(dp),    intent(in)     :: T_adim,Tmax_adim
    integer(sp), intent(inout)  :: aux_matrix_pbc(n+2,n+2)
    real(dp),    intent(inout)  :: U_adim       ! Energía interna
    real(dp),    intent(inout)  :: Madim ! Magenitación
    real(dp),    intent(inout)  :: autocor_vector_mom_U(tau_corr),autocor_vector_mom_M(tau_corr)

    integer(sp), parameter   :: m_exp=3_sp ! numero de experimentos
    integer(sp)              :: i,j
    real(dp),    allocatable :: autocor_vector_U(:),aux_vector1_U(:),aux_vector2_U(:),mask_vector_U(:)
    real(dp),    allocatable :: autocor_vector_M(:),aux_vector1_M(:),aux_vector2_M(:),mask_vector_M(:)
    real(dp)                 :: U_med,var_U,M_med,var_M ! valor media y varianza (autocorrelación)

    allocate(autocor_vector_U(tau_corr),aux_vector1_U(tau_corr),&
             aux_vector2_U(tau_corr),mask_vector_U(tau_corr))
    allocate(autocor_vector_M(tau_corr),aux_vector1_M(tau_corr),&
             aux_vector2_M(tau_corr),mask_vector_M(tau_corr))
    if (T_adim==Tmax_adim) then
        autocor_vector_mom_U(:)=0._dp;U_med=0._dp;var_U=0._dp
        autocor_vector_mom_M(:)=0._dp;M_med=0._dp;var_M=0._dp
    end if
    do j=1,m_exp
        do i=1,MC_step
            call MC_step_relaxation(1_sp,n,aux_matrix_pbc,T_adim,U_adim)
            Madim=M_adim(aux_matrix_pbc,n) ! Magnetización
            ! datos para hacer estadística en steady state
            if (i>=MC_step_trans.and.T_adim==Tmax_adim) then! calcular autocorr para la Temp requerida             
                call func_autocor(U_adim,i-MC_step_trans+1,tau_corr,num_of_terms,&
                                autocor_vector_U,aux_vector1_U,mask_vector_U,aux_vector2_U,&
                                U_med,var_U)
                call func_autocor(Madim,i-MC_step_trans+1,tau_corr,num_of_terms,&
                                autocor_vector_M,aux_vector1_M,mask_vector_M,aux_vector2_M,&
                                M_med,var_M)
            end if
        end do
        !autocorrelación
        if (T_adim==Tmax_adim) then
            autocor_vector_U(:)=(autocor_vector_U(:)-U_med*U_med)*(1._dp/var_U)
            autocor_vector_mom_U(:)=autocor_vector_mom_U(:)+autocor_vector_U(:)
            autocor_vector_M(:)=(autocor_vector_M(:)-M_med*M_med)*(1._dp/var_M)
            autocor_vector_mom_M(:)=autocor_vector_mom_M(:)+autocor_vector_M(:)
        end if
    end do
    deallocate(aux_vector1_U,aux_vector2_U,mask_vector_U)
    deallocate(aux_vector1_M,aux_vector2_M,mask_vector_M)
    if (T_adim==Tmax_adim) then
        autocor_vector_mom_U(:)=(1._dp/real(m_exp,dp))*autocor_vector_mom_U(:)
        autocor_vector_mom_M(:)=(1._dp/real(m_exp,dp))*autocor_vector_mom_M(:)
    end if
    deallocate(autocor_vector_U,autocor_vector_M)
end subroutine rlx_ising