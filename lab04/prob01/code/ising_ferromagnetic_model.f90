! P01.a
! make clean && make ising_ferromagnetic_model.o && ./ising_ferromagnetic_model.o
program ising_ferromagnetic_model
    use module_precision;use module_2D_ferromagnetic_ising
    implicit none

    integer(sp), parameter   :: n=40_sp,m=100    ! sitios de red por dimension
    integer(sp), parameter   :: MC_step=50000_sp,MC_step_trans=3000_sp  ! Monte Carlo step total and transitory
    real(dp),    parameter   :: Tc_adim=2.2676_dp ! temperatura de Curie adimensional
    integer(sp), allocatable :: aux_matrix_pbc(:,:)
    real(dp)                 :: T_adim,T_step ! temperatura adimensional
    integer(sp)              :: i,j,istat
    real(dp)                 :: U_adim,U_med_adim,sigma_U,error_U
    real(dp)                 :: Madim,M_med_adim,sigma_M,error_M
    real(dp)                 :: s0,s1_U,s2_U,s1_M,s2_M ! variables para hacer estadística

    open(10,file='../results/result_01a_MCS_trans_T0.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat
    21 format(2(A12,x),A12);write(10,21) 'MC_step','U_adim','M_adim'
    20 format(I12,x,E12.4,x,E12.4)

    allocate(aux_matrix_pbc(n+2,n+2))
    ! genero configuracion inicial random
    call initial_spins_configuration(1_sp,n,aux_matrix_pbc)
    ! calculamos energía interna (configuración inicial)
    call average_energy(aux_matrix_pbc,n,U_adim)
    ! calculamos magnetización inicial
    Madim=M_adim(aux_matrix_pbc,n)

    ! temperaturas de equilibrio
    !T0=1;T1=2.2;T2=2.2676;T3=3.3

    s0=0._dp;s1_U=0._dp;s2_U=0._dp;s1_M=0._dp;s2_M=0._dp

    do1: do j=1,m ! cantidad de Temperaturas diferentes
        T_step=3.3_dp*(1._dp/real(m-1_sp,dp))
        T_adim=T_step*real(j-1_sp,dp)
        write(*,'(A12,E12.4)') 'T_adim=',T_adim
        if (T_adim==1._dp) write(10,20) 0,U_adim,Madim
        do i=1,MC_step
            call MC_step_relaxation(1_sp,n,aux_matrix_pbc,T_adim,U_adim)
            Madim=M_adim(aux_matrix_pbc,n) ! Magnetización
            if (T_adim==1._dp) write(10,20) i,U_adim,Madim

            ! datos para hacer estadística en steady state
            ! if (i>=MC_step_trans) then
            !     s0=s0+1._dp
            !     s2_U=s2_U+U_adim*U_adim;s1_U=s1_U+U_adim ! energía
            !     s2_M=s2_M+Madim*Madim;s1_M=s1_M+Madim    ! magnetización
            ! end if
        end do
        if (T_adim==1._dp) exit do1
    end do do1
    close(10)

    ! !calculamos media,desviacion estándar, varianza y error
    ! U_med_adim=s1_U*(1._dp/s0)
    ! sigma_U=sqrt((s2_U-s0*U_med_adim*U_med_adim)*(1._dp/(s0-1._dp)))
    ! error_U=sigma_U*(1._dp/sqrt(s0-1._dp))
    ! M_med_adim=s1_M*(1._dp/s0)
    ! sigma_M=sqrt((s2_M-s0*M_med_adim*M_med_adim)*(1._dp/(s0-1._dp)))
    ! error_M=sigma_M*(1._dp/sqrt(s0-1._dp))

    ! write(*,'(A14,E12.4)') 'M_exact_adim=', M_exact_adim(n,T_adim(3),Tc_adim)
    ! write(*,'(2(A14,E12.4))') 'M_med_adim=', M_med_adim,'error_M',error_M
    ! write(*,'(2(A14,E12.4))') 'U_med_adim=', U_med_adim,'error_U',error_U

    deallocate(aux_matrix_pbc)
end program ising_ferromagnetic_model