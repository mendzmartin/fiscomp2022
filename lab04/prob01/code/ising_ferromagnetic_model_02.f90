! P01.b (TyMEII G04.P16;G04.P02)
! make clean && make ising_ferromagnetic_model_02.o && ./ising_ferromagnetic_model_02.o
program ising_ferromagnetic_model_02
    use module_precision;use module_2D_ferromagnetic_ising
    implicit none

    integer(sp), parameter   :: n=128_sp                                 ! sitios de red por dimension
    integer(sp), parameter   :: MC_step=10000_sp,MC_step_trans=3000_sp  ! Monte Carlo step total and transitory
    integer(sp), parameter   :: m1=100_sp,m2=1000_sp                    ! puntos p/ mallado fino y grueso
    real(dp),    parameter   :: Tmin_adim=0._dp,Tmax_adim=10._dp!7.5_dp        ! temperatura adimensional
    real(dp),    parameter   :: Tc_adim=2.2676_dp                       ! temperatura de Curie adimensional
    real(dp),    parameter   :: deltaT_adim=2.2676_dp*0.5_dp            ! intervalo para incremetar ptos
    integer(sp), allocatable :: aux_matrix_pbc(:,:)
    integer(sp)              :: i,istat
    real(dp)                 :: U_adim,U_med_adim,sigma_U,error_U       ! Energía interna
    real(dp)                 :: Madim,M_med_adim,sigma_M,error_M,Mexact ! Magenitación
    real(dp)                 :: susc_adim                               ! Susceptibilidad magnética
    real(dp)                 :: cv                                      ! calor específico 
    real(dp)                 :: T_adim,T_step1,T_step2

    open(10,file='../results/result_01b_07.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat
    21 format(7(A12,x),A12);write(10,21) 'T/Tc','U_med','U_error','M_med','M_error','Mexact','cv','susc'
    20 format(7(E12.4,x),E12.4)

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
                       susc_adim,cv)
        write(10,20) T_adim*(1._dp/Tc_adim),U_med_adim,error_U,M_med_adim,error_M,Mexact,cv,susc_adim
        write(*,*) i
        !if (T_adim==0._dp) write(*,*) 'U_med_adim=',U_med_adim
    end do

    ! mallado fino
    do i=1,m2
        T_step2=2_dp*deltaT_adim*(1._dp/real(m2-1,dp))
        T_adim=(Tc_adim-deltaT_adim)+T_step2*real(i-1,dp)
        call rlx_ising(n,aux_matrix_pbc,MC_step,MC_step_trans,T_adim,Tc_adim,&
                       U_adim,U_med_adim,sigma_U,error_U,&
                       Madim,M_med_adim,sigma_M,error_M,Mexact,&
                       susc_adim,cv)
        write(10,20) T_adim*(1._dp/Tc_adim),U_med_adim,error_U,M_med_adim,error_M,Mexact,cv,susc_adim
        write(*,*) i+m1
    end do

    ! mallado grueso
    do i=1,m1
        T_step1=abs(Tmax_adim-(Tc_adim+deltaT_adim))*(1._dp/real(m1-1,dp))
        T_adim=(Tc_adim+deltaT_adim)+T_step1*real(i-1,dp)
        call rlx_ising(n,aux_matrix_pbc,MC_step,MC_step_trans,T_adim,Tc_adim,&
                       U_adim,U_med_adim,sigma_U,error_U,&
                       Madim,M_med_adim,sigma_M,error_M,Mexact,&
                       susc_adim,cv)
        write(10,20) T_adim*(1._dp/Tc_adim),U_med_adim,error_U,M_med_adim,error_M,Mexact,cv,susc_adim
        write(*,*) i+m1+m2
    end do

    !close(10)
    deallocate(aux_matrix_pbc)
end program ising_ferromagnetic_model_02

subroutine rlx_ising(n,aux_matrix_pbc,MC_step,MC_step_trans,T_adim,Tc_adim,&
                     U_adim,U_med_adim,sigma_U,error_U,&
                     Madim,M_med_adim,sigma_M,error_M,Mexact,&
                     susc_adim,cv)
    use module_precision;use module_2D_ferromagnetic_ising
    implicit none

    integer(sp), intent(in)     :: n,MC_step,MC_step_trans
    real(dp),    intent(in)     :: T_adim,Tc_adim
    integer(sp), intent(inout)  :: aux_matrix_pbc(n+2,n+2)
    real(dp),    intent(inout)  :: U_adim,U_med_adim,sigma_U,error_U       ! Energía interna
    real(dp),    intent(inout)  :: Madim,M_med_adim,sigma_M,error_M,Mexact        ! Magenitación
    real(dp),    intent(inout)  :: susc_adim                               ! Susceptibilidad magnética
    real(dp),    intent(inout)  :: cv                                      ! calor específico 

    real(dp)    :: s0,s1_U,s2_U,s1_M,s2_M ! variables para hacer estadística
    integer(sp) :: i

    s0=0._dp;s1_U=0._dp;s2_U=0._dp;s1_M=0._dp;s2_M=0._dp
    do i=1,MC_step
        call MC_step_relaxation(1_sp,n,aux_matrix_pbc,T_adim,U_adim)
        Madim=M_adim(aux_matrix_pbc,n) ! Magnetización
        ! datos para hacer estadística en steady state
        if (i>=MC_step_trans) then
            s0=s0+1._dp
            s2_U=s2_U+U_adim*U_adim;s1_U=s1_U+U_adim ! energía
            s2_M=s2_M+Madim*Madim;s1_M=s1_M+Madim    ! magnetización
        end if
    end do

    ! calculamos media,desviacion estándar, varianza y error
    U_med_adim=s1_U*(1._dp/s0)
    sigma_U=sqrt((s2_U-s0*U_med_adim*U_med_adim)*(1._dp/(s0-1._dp)))
    error_U=sigma_U*(1._dp/sqrt(s0-1._dp))
    M_med_adim=s1_M*(1._dp/s0)
    sigma_M=sqrt((s2_M-s0*M_med_adim*M_med_adim)*(1._dp/(s0-1._dp)))
    error_M=sigma_M*(1._dp/sqrt(s0-1._dp))
    if (T_adim==0._dp) then;cv=0._dp;susc_adim=0._dp
    else
        cv=(1._dp/T_adim)*sigma_U*sigma_U
        susc_adim=(1._dp/T_adim)*sigma_M*sigma_M
    end if
    Mexact=M_exact_adim(n,T_adim,Tc_adim)

end subroutine rlx_ising