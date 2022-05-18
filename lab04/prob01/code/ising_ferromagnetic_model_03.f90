! Programa para determinar el tiempo de termalización
!  se realizan corridas del método de MC "metropolis" para distintas temperaturas
!  tanto cerca como lejos de la temperatura crítica.
! make clean && make ising_ferromagnetic_model_03.o && ./ising_ferromagnetic_model_03.o
program ising_ferromagnetic_model_03
    use module_precision;use module_2D_ferromagnetic_ising
    implicit none

    integer(sp), parameter   :: n=40_sp,m=10_sp     ! sitios de red por dimension
    integer(sp), parameter   :: MC_step=10000_sp    ! Monte Carlo step total and transitory
    real(dp),    parameter   :: Tc_adim=2.2676_dp   ! temperatura de Curie adimensional
    integer(sp), allocatable :: aux_matrix_pbc(:,:)
    integer(sp)              :: istat
    real(dp)                 :: Madim,U_adim

    allocate(aux_matrix_pbc(n+2,n+2))
    ! genero configuracion inicial random
    call initial_spins_configuration(1_sp,n,aux_matrix_pbc)
    ! calculamos energía interna (configuración inicial)
    call average_energy(aux_matrix_pbc,n,U_adim)
    ! calculamos magnetización inicial
    Madim=M_adim(aux_matrix_pbc,n)

    ! temperaturas de equilibrio
    !T0=1;T1=2.2;T2=2.2676;T3=3.3

    21 format(2(A12,x),A12)

    open(10,file='../results/result_01a_MCS_trans_T0.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat
    write(10,21) 'MC_step','U_adim','M_adim'
    call ising_relax(n,MC_step,m,10,aux_matrix_pbc,0._dp,1._dp,U_adim,Madim)
    close(10)

    open(11,file='../results/result_01a_MCS_trans_T1.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(11file) = ',istat
    write(11,21) 'MC_step','U_adim','M_adim'
    call ising_relax(n,MC_step,m,11,aux_matrix_pbc,1.1_dp,2._dp,U_adim,Madim)
    close(11)

    open(12,file='../results/result_01a_MCS_trans_T2.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(12file) = ',istat
    write(12,21) 'MC_step','U_adim','M_adim'
    call ising_relax(n,MC_step,m,12,aux_matrix_pbc,2.1_dp,Tc_adim,U_adim,Madim)
    close(12)

    open(13,file='../results/result_01a_MCS_trans_T3.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(13file) = ',istat
    write(13,21) 'MC_step','U_adim','M_adim'
    call ising_relax(n,MC_step,m,13,aux_matrix_pbc,Tc_adim+0.1_dp,3.3_dp,U_adim,Madim)
    close(13)

    deallocate(aux_matrix_pbc)
end program ising_ferromagnetic_model_03

subroutine ising_relax(n,MC_step,m,file_num,aux_matrix_pbc,T_start,T_end,U_adim,Madim)
    use module_precision;use module_2D_ferromagnetic_ising
    implicit none 
    real(dp),    intent(in)     :: T_end,T_start
    integer(sp), intent(in)     :: n,m,file_num,MC_step
    integer(sp), intent(inout)  :: aux_matrix_pbc(n+2,n+2)
    real(dp),    intent(inout)  :: U_adim,Madim
    integer(sp)                 :: i,j
    real(dp)                    :: T_step,T_adim

    20 format(I12,x,E12.4,x,E12.4)

    bucle_01: do j=1,m ! cantidad de Temperaturas diferentes
        T_step=abs(T_end-T_start)*(1._dp/real(m-1_sp,dp))
        T_adim=T_start+T_step*real(j-1_sp,dp)
        write(*,'(A12,E12.4)') 'T_adim=',T_adim
        do i=1,MC_step
            call MC_step_relaxation(1_sp,n,aux_matrix_pbc,T_adim,U_adim)
            Madim=M_adim(aux_matrix_pbc,n) ! Magnetización
            if (T_adim==T_end) write(file_num,20) i,U_adim,Madim
        end do
    end do bucle_01

end subroutine ising_relax