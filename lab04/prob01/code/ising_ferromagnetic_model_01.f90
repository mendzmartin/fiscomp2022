! Programa para resolver el inciso a)
! (primero se uso el ising_ferromagnetic_model_03.f90 para determinar MC_step_trans)
! make clean && make ising_ferromagnetic_model_01.o && ./ising_ferromagnetic_model_01.o
program ising_ferromagnetic_model_01
    use module_precision;use module_2D_ferromagnetic_ising
    implicit none

    integer(sp), parameter   :: n=40_sp,m=10_sp                        ! sitios de red por dimension
    integer(sp), parameter   :: MC_step=10000_sp,MC_step_trans=3000_sp ! Monte Carlo step total and transitory
    real(dp),    parameter   :: Tc_adim=2.2676_dp                      ! temperatura de Curie adimensional
    integer(sp), allocatable :: aux_matrix_pbc(:,:)
    integer(sp)              :: i,j,istat
    real(dp)                 :: Madim,U_adim

    allocate(aux_matrix_pbc(n+2,n+2))
    
    !call initial_spins_configuration(2_sp,n,aux_matrix_pbc) ! genero configuracion inicial random
    call initial_spins_configuration(1_sp,n,aux_matrix_pbc) ! genero configuracion inicial todos up

    ! calculamos energía interna (configuración inicial)
    call average_energy(aux_matrix_pbc,n,U_adim)
    ! calculamos magnetización inicial
    Madim=M_adim(aux_matrix_pbc,n)

    ! temperaturas de equilibrio
    !T0=1;T1=2.2;T2=2.2676;T3=3.3

    21 format(2(A12,x),A12)

    open(10,file='../results/result_01a_T=1_up.dat',status='replace',action='write',iostat=istat)  ! ordenado
    !open(10,file='../results/result_01a_T=1_rnd.dat',status='replace',action='write',iostat=istat) ! desordenado
    write(*,*) 'istat(10file) = ',istat
    write(10,21) 'MC_step','U_adim','M_adim'
    call ising_relax(n,MC_step,MC_step_trans,m,10,aux_matrix_pbc,0._dp,1._dp,Tc_adim,U_adim,Madim)
    close(10)

    ! mostramos mapa de spins antes de Tc
    open(50,file='../results/result_01a_spinmap_T=1.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(50file) = ',istat
    do j=2,n+1;do i=2,n+1
        write(50,'(2(I5,x),I5)') i-1,j-1,aux_matrix_pbc(i,j)
    end do;end do
    close(50)

    open(11,file='../results/result_01a_T=2_up.dat',status='replace',action='write',iostat=istat)  ! ordenado
    !open(11,file='../results/result_01a_T=2_rnd.dat',status='replace',action='write',iostat=istat) ! desordenado
    write(*,*) 'istat(11file) = ',istat
    write(11,21) 'MC_step','U_adim','M_adim'
    call ising_relax(n,MC_step,MC_step_trans,m,11,aux_matrix_pbc,1.1_dp,2._dp,Tc_adim,U_adim,Madim)
    close(11)

    open(12,file='../results/result_01a_T=Tc_up.dat',status='replace',action='write',iostat=istat)  ! ordenado
    !open(12,file='../results/result_01a_T=Tc_rnd.dat',status='replace',action='write',iostat=istat) ! desordenado
    write(*,*) 'istat(12file) = ',istat
    write(12,21) 'MC_step','U_adim','M_adim'
    call ising_relax(n,MC_step,MC_step_trans,m,12,aux_matrix_pbc,2.1_dp,Tc_adim,Tc_adim,U_adim,Madim)
    close(12)

    ! mostramos mapa de spins en Tc
    open(50,file='../results/result_01a_spinmap_T=Tc.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(50file) = ',istat
    do j=2,n+1;do i=2,n+1
        write(50,'(2(I5,x),I5)') i-1,j-1,aux_matrix_pbc(i,j)
    end do;end do
    close(50)


    open(13,file='../results/result_01a_T=3.3_up.dat',status='replace',action='write',iostat=istat)  ! ordenado
    !open(13,file='../results/result_01a_T=3.3_rnd.dat',status='replace',action='write',iostat=istat) ! desordenado
    write(*,*) 'istat(13file) = ',istat
    write(13,21) 'MC_step','U_adim','M_adim'
    call ising_relax(n,MC_step,MC_step_trans,m,13,aux_matrix_pbc,Tc_adim+0.1_dp,3.3_dp,Tc_adim,U_adim,Madim)
    close(13)

    ! mostramos mapa de spins luego de Tc
    open(50,file='../results/result_01a_spinmap_T=3.3.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(50file) = ',istat
    do j=2,n+1;do i=2,n+1
        write(50,'(2(I5,x),I5)') i-1,j-1,aux_matrix_pbc(i,j)
    end do;end do
    close(50)

    deallocate(aux_matrix_pbc)
end program ising_ferromagnetic_model_01

subroutine ising_relax(n,MC_step,MC_step_trans,m,file_num,aux_matrix_pbc,T_start,T_end,Tc_adim,U_adim,Madim)
    use module_precision;use module_2D_ferromagnetic_ising
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
            call MC_step_relaxation(1_sp,n,aux_matrix_pbc,T_adim,U_adim)
        end do
    end do
    ! termalizo a la temperatura buscada
    write(*,'(A12,E12.4)') 'T_adim=',T_end
    s0=0._dp;s1_U=0._dp;s2_U=0._dp;s1_M=0._dp;s2_M=0._dp
    do i=1,MC_step
        call MC_step_relaxation(1_sp,n,aux_matrix_pbc,T_end,U_adim)
        Madim=M_adim(aux_matrix_pbc,n) ! Magnetización
        ! datos para hacer estadística en steady state
        if (i>=MC_step_trans) then
            s0=s0+1._dp
            s2_U=s2_U+U_adim*U_adim;s1_U=s1_U+U_adim ! energía
            s2_M=s2_M+Madim*Madim;s1_M=s1_M+Madim    ! magnetización
        end if
        write(file_num,20) i,U_adim,Madim
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

end subroutine ising_relax