! Programa para calcular histogramas (inciso c)
! make clean && make ising_ferromagnetic_model_04.o && ./ising_ferromagnetic_model_04.o
program ising_ferromagnetic_model_04
    use module_precision;use module_2D_ferromagnetic_ising
    implicit none

    integer(sp), parameter   :: MC_step=10000_sp,MC_step_rlx=3000_sp     ! Monte Carlo step total and transitory
    real(dp),    parameter   :: Tc_adim=2.2692_dp    ! temperatura de Curie adimensional
    integer(sp), allocatable :: aux_matrix_pbc(:,:)
    real(dp),    allocatable :: Madim_vector(:)    ! vector de magnetización para histograma
    integer(sp)              :: n                    ! sitios de red por dimension
    integer(sp)              :: m              ! cantidad puntos para definir el paso de temperaturas
    real(dp)                 :: Madim,U_adim


    m=10_sp;n=10_sp;allocate(aux_matrix_pbc(n+2,n+2));allocate(Madim_vector(MC_step-MC_step_rlx))
    call initial_spins_configuration(1_sp,n,aux_matrix_pbc)
    call average_energy(aux_matrix_pbc,n,U_adim);Madim=M_adim(aux_matrix_pbc,n)
    call ising_relax(n,MC_step,MC_step_rlx,m,aux_matrix_pbc,0._dp,Tc_adim*1.1999_dp,U_adim,Madim,Madim_vector)
    call histogram('../results/result_01c_historgram_lessTc_10x10.dat',10_sp,Madim_vector,MC_step-MC_step_rlx)
    call initial_spins_configuration(1_sp,n,aux_matrix_pbc)
    call average_energy(aux_matrix_pbc,n,U_adim);Madim=M_adim(aux_matrix_pbc,n)
    call ising_relax(n,MC_step,MC_step_rlx,m,aux_matrix_pbc,0._dp,Tc_adim*1.5_dp,U_adim,Madim,Madim_vector)
    call histogram('../results/result_01c_historgram_greaterTc_10x10.dat',10_sp,Madim_vector,MC_step-MC_step_rlx)
    deallocate(aux_matrix_pbc);deallocate(Madim_vector)

    m=10_sp;n=20_sp;allocate(aux_matrix_pbc(n+2,n+2));allocate(Madim_vector(MC_step-MC_step_rlx))
    call initial_spins_configuration(1_sp,n,aux_matrix_pbc)
    call average_energy(aux_matrix_pbc,n,U_adim);Madim=M_adim(aux_matrix_pbc,n)
    call ising_relax(n,MC_step,MC_step_rlx,m,aux_matrix_pbc,0._dp,Tc_adim*1.11_dp,U_adim,Madim,Madim_vector)
    call histogram('../results/result_01c_historgram_lessTc_20x20.dat',10_sp,Madim_vector,MC_step-MC_step_rlx)
    call initial_spins_configuration(1_sp,n,aux_matrix_pbc)
    call average_energy(aux_matrix_pbc,n,U_adim);Madim=M_adim(aux_matrix_pbc,n)
    call ising_relax(n,MC_step,MC_step_rlx,m,aux_matrix_pbc,0._dp,Tc_adim*1.5_dp,U_adim,Madim,Madim_vector)
    call histogram('../results/result_01c_historgram_greaterTc_20x20.dat',10_sp,Madim_vector,MC_step-MC_step_rlx)
    deallocate(aux_matrix_pbc);deallocate(Madim_vector)

    m=10_sp;n=40_sp;allocate(aux_matrix_pbc(n+2,n+2));allocate(Madim_vector(MC_step-MC_step_rlx))
    call initial_spins_configuration(1_sp,n,aux_matrix_pbc)
    call average_energy(aux_matrix_pbc,n,U_adim);Madim=M_adim(aux_matrix_pbc,n)
    call ising_relax(n,MC_step,MC_step_rlx,m,aux_matrix_pbc,0._dp,Tc_adim*1.118_dp,U_adim,Madim,Madim_vector)
    call histogram('../results/result_01c_historgram_lessTc_40x40.dat',10_sp,Madim_vector,MC_step-MC_step_rlx)
    call initial_spins_configuration(1_sp,n,aux_matrix_pbc)
    call average_energy(aux_matrix_pbc,n,U_adim);Madim=M_adim(aux_matrix_pbc,n)
    call ising_relax(n,MC_step,MC_step_rlx,m,aux_matrix_pbc,0._dp,Tc_adim*1.5_dp,U_adim,Madim,Madim_vector)
    call histogram('../results/result_01c_historgram_greaterTc_40x40.dat',10_sp,Madim_vector,MC_step-MC_step_rlx)
    deallocate(aux_matrix_pbc);deallocate(Madim_vector)

end program ising_ferromagnetic_model_04

subroutine ising_relax(n,MC_step,MC_step_rlx,m,aux_matrix_pbc,T_start,T_end,U_adim,Madim,Madim_vector)
    use module_precision;use module_2D_ferromagnetic_ising
    implicit none 
    real(dp),    intent(in)     :: T_end,T_start
    integer(sp), intent(in)     :: n,m,MC_step,MC_step_rlx
    integer(sp), intent(inout)  :: aux_matrix_pbc(n+2,n+2)
    real(dp),    intent(inout)  :: U_adim,Madim
    real(dp),    intent(inout)    :: Madim_vector(MC_step-MC_step_rlx)
    integer(sp)                 :: i,j
    real(dp)                    :: T_step,T_adim

    bucle_01: do j=1,m ! cantidad de Temperaturas diferentes
        T_step=abs(T_end-T_start)*(1._dp/real(m-1_sp,dp))
        T_adim=T_start+T_step*real(j-1_sp,dp)
        do i=1,MC_step_rlx
            call MC_step_relaxation(1_sp,n,aux_matrix_pbc,T_adim,U_adim)
            Madim=M_adim(aux_matrix_pbc,n) ! Magnetización
        end do
        do i=MC_step_rlx+1,MC_step
            call MC_step_relaxation(1_sp,n,aux_matrix_pbc,T_adim,U_adim)
            Madim=M_adim(aux_matrix_pbc,n) ! Magnetización
            if (j==m) Madim_vector(i-MC_step_rlx)=Madim
        end do
    end do bucle_01

end subroutine ising_relax

! subrutina para crear e imprimir histograma
subroutine histogram(file_name,file_number,x_vector,x_dim)
    use module_precision

    implicit none
    character(len=*), intent(in) :: file_name
    integer(sp),      intent(in) :: file_number
    integer(sp),      intent(in) :: x_dim           ! dimension
    real(dp),         intent(in) :: x_vector(x_dim) ! normalized data

    integer(sp), parameter   :: n_bins=1000_sp   ! numbers of bins (JUGAR CON ESTE VALOR)
    real(dp),    allocatable :: bins_points(:)  ! points between bins vector
    integer(sp), allocatable :: counter(:)      ! counter vector of bins
    real(sp)                 :: max_value       ! maximun counter value
    real(dp)                 :: min_bin_point,max_bin_point
    real(dp)                 :: bins_step       ! step of points between bins
    integer(sp)              :: i,j,istat       ! loop and control variables

    ! armamos el vector de bins
    allocate(bins_points(n_bins+1),counter(n_bins))
    max_bin_point=1._dp;min_bin_point=-1_dp
    bins_step=abs(max_bin_point-min_bin_point)*(1._dp/n_bins)
    do i=1,n_bins+1;bins_points(i)=min_bin_point+bins_step*real(i-1,dp);end do

    ! llenamos el vector contador de bins
    do i=1,n_bins
        counter(i)=0
        do j=1,x_dim
            if ((bins_points(i)<=x_vector(j)).and.(bins_points(i+1)>=x_vector(j))) then
                counter(i)=counter(i)+1
            endif
    enddo;enddo

    max_value=real(maxval(counter(:)),dp)
    if (max_value==0._dp) write(*,*) 'math error'

    ! escribimos el histograma
    open(file_number,file=file_name,status='replace',action='write',iostat=istat)
    21 format(A12,x,A12); 20 format(E12.4,x,E12.4);write(*,*) 'istat=', istat
    write(file_number,21) 'bins points','counter'
    do i = 1,n_bins; write(file_number,20) bins_points(i),real(counter(i),dp)*(1._dp/max_value); enddo
    close(file_number)
end subroutine histogram