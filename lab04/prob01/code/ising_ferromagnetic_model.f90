! P01.a
program ising_ferromagnetic_model
    use module_precision;use module_mt19937, only: sgrnd,grnd
    implicit none

    integer(sp), parameter   :: n=10_sp           ! sitios de red por dimension
    integer(sp), parameter   :: MC_step=10000_sp  ! Monte Carlo step
    integer(sp), parameter   :: T_init_type=1_sp,MC_step_type=1._sp
    real(dp),    parameter   :: Tc_adim=2.2676_dp ! temperatura de Curie adimensional
    integer(sp), allocatable :: spin_matrix(:,:),aux_matrix_pbc(:,:)
    real(dp)                 :: T_adim ! temperatura adimensional
    integer(sp)              :: seed,seed_val(8),i,j,k,l,istat
    real(dp)                 :: nrand
    real(dp)                 :: U_adim,s0,s1_U,s2_U,sigma_U,var_U,error_U
    integer(sp)              :: U_delta_adim,deltaU_adim
    real(dp)                 :: Madim,M_adim,M_exact_adim,s1_M,s2_M,sigma_M,var_M,error_M

    open(10,file='../results/result_01a.dat',status='replace',action='write',iostat=istat)
    write(*,*) 'istat(10file) = ',istat
    21 format(4(A12,x),A12); write(10,21) 'MC_step','U_adim','U_error','M_adim','M_error'
    20 format(I12,x,3(E12.4,x),E12.4)

    call date_and_time(values=seed_val)
    seed=seed_val(8)*seed_val(7)*seed_val(6)+seed_val(5);call sgrnd(seed)

    allocate(spin_matrix(n,n))
    
    select case (T_init_type)
        case(1) ! genero configuracion inicial (random ó T(inicial)=Infinity)
            do j=1,n;do i=1,n
                nrand=real(grnd(),dp)
                if (nrand>=0._dp.and.nrand<0.5_dp) then;spin_matrix(i,j)=1_sp
                else; spin_matrix(i,j)=-1_sp;end if
            end do;end do
        case(2) ! genero configuracion inicial (ordenada ó T(inicial)=0)
            spin_matrix(:,:) = 1_sp ! todos los spins "up"
    end select

    allocate(aux_matrix_pbc(n+2,n+2))
    ! creamos matriz de spin con PBC
    aux_matrix_pbc(1,2:n+1)=spin_matrix(n,1:n)
    aux_matrix_pbc(n+2,2:n+1)=spin_matrix(1,1:n)
    aux_matrix_pbc(2:n+1,1)=spin_matrix(1:n,n)
    aux_matrix_pbc(2:n+1,n+2)=spin_matrix(1:n,1)
    aux_matrix_pbc(2:n+1,2:n+1)=spin_matrix(1:n,1:n)
    deallocate(spin_matrix)

    ! calculamos energía interna (configuración inicial)
    call average_energy(aux_matrix_pbc,n,U_adim)
    write(10,20) 0_sp,U_adim,0._dp,M_adim(aux_matrix_pbc,n),0._dp

    T_adim=2.26_dp ! T_adim={2.00;3.30;2.26}
    call date_and_time(values=seed_val)
    seed=seed_val(8)*seed_val(7)*seed_val(6)+seed_val(5);call sgrnd(seed)

    s0=0._dp;s1_U=0._dp;s2_U=0._dp
    s1_M=0._dp;s2_M=0._dp

    select case(MC_step_type)
        case(1) ! hago un MC_step de forma random
        do k=1,MC_step
            do l=1,n*n
                ! floor(x,type) returns the greatest integer less than or equal to X
                nrand=real(grnd(),dp);i=2_sp+floor(n*nrand,sp) ! random integer from 2 to n+2
                nrand=real(grnd(),dp);j=2_sp+floor(n*nrand,sp)
                deltaU_adim=U_delta_adim(aux_matrix_pbc(i,j),&
                                        aux_matrix_pbc(i,j+1),aux_matrix_pbc(i,j-1),&
                                        aux_matrix_pbc(i+1,j),aux_matrix_pbc(i-1,j))
                if (deltaU_adim<0_sp) then
                    aux_matrix_pbc(i,j)=-aux_matrix_pbc(i,j)
                    U_adim=U_adim+deltaU_adim
                else 
                    nrand=real(grnd(),dp)
                    if (exp(-deltaU_adim*(1._dp/T_adim))>=nrand) then
                        aux_matrix_pbc(i,j)=-aux_matrix_pbc(i,j)
                        U_adim=U_adim+deltaU_adim
                end if;end if
            end do
            Madim=M_adim(aux_matrix_pbc,n)
            ! calculamos desviacion estándar, varianza y error
            s0=s0+real(n*n,dp)
            ! energía
            s2_U=s2_U+U_adim*U_adim;s1_U=s1_U+U_adim
            sigma_U=sqrt(s0*s2_U-(s1_U)*(s1_U))*(1._dp/s0);var_U=sigma_U*sigma_U
            error_U=sigma_U*(1._dp/sqrt(s0))
            ! magnetización
            s2_M=s2_M+Madim*Madim;s1_M=s1_M+Madim
            sigma_M=sqrt(s0*s2_M-(s1_M)*(s1_M))*(1._dp/s0);var_M=sigma_M*sigma_M
            error_M=sigma_M*(1._dp/sqrt(s0))
            write(10,20) k,U_adim,error_U,Madim,error_M
        end do
        case(2) ! hago un paso de Monte Carlo (MC_step) de forma ordenada
        do k=1,MC_step
            do j=2,n+1;do i=2,n+1
                deltaU_adim=U_delta_adim(aux_matrix_pbc(i,j),&
                                        aux_matrix_pbc(i,j+1),aux_matrix_pbc(i,j-1),&
                                        aux_matrix_pbc(i+1,j),aux_matrix_pbc(i-1,j))
                if (deltaU_adim<0_sp) then
                    aux_matrix_pbc(i,j)=-aux_matrix_pbc(i,j)
                    U_adim=U_adim+deltaU_adim
                else 
                    nrand=real(grnd(),dp)
                    if (exp(-deltaU_adim*(1._dp/T_adim))>=nrand) then
                        aux_matrix_pbc(i,j)=-aux_matrix_pbc(i,j)
                        U_adim=U_adim+deltaU_adim
                end if;end if
            end do;end do
            Madim=M_adim(aux_matrix_pbc,n)
            ! calculamos desviacion estándar, varianza y error
            s0=s0+real(n*n,dp)
            ! energía
            s2_U=s2_U+U_adim*U_adim;s1_U=s1_U+U_adim
            sigma_U=sqrt(s0*s2_U-(s1_U)*(s1_U))*(1._dp/s0);var_U=sigma_U*sigma_U
            error_U=sigma_U*(1._dp/sqrt(s0))
            ! magnetización
            s2_M=s2_M+Madim*Madim;s1_M=s1_M+Madim
            sigma_M=sqrt(s0*s2_M-(s1_M)*(s1_M))*(1._dp/s0);var_M=sigma_M*sigma_M
            error_M=sigma_M*(1._dp/sqrt(s0))
            write(10,20) k,U_adim,error_U,Madim,error_M
        end do
    end select
    close(10)

    write(*,'(A14,E12.4)') 'M_exact_adim=',M_exact_adim(n,T_adim,Tc_adim)

    deallocate(aux_matrix_pbc)
end program ising_ferromagnetic_model

! Energia en unidades de J_int (coeficiente de interacción >0 => ferromagnetic)
subroutine average_energy(aux_matrix_pbc,n,U)
    use module_precision
    implicit none
    integer(sp), intent(in)    :: n
    integer(sp), intent(in)    :: aux_matrix_pbc(1:n+2,1:n+2)
    real(dp),    intent(inout) :: U
    integer(sp) :: i,j

    U=0._dp
    do j=2,n
        do i=2,n
            U=U+real(aux_matrix_pbc(i,j)*(aux_matrix_pbc(i,j+1)+aux_matrix_pbc(i,j-1)&
                   & +aux_matrix_pbc(i+1,j)+aux_matrix_pbc(i-1,j)),dp)
        end do
    end do
    U=-1._dp*U ! suponemos J_int >0
end subroutine average_energy

! Delta energético por "flipear" un spin, en unidades de J_int
function U_delta_adim(spin_value,spin1,spin2,spin3,spin4)
    use module_precision
    implicit none
    integer(sp), intent(in) :: spin_value,spin1,spin2,spin3,spin4
    integer(sp)             :: val,U_delta_adim
    val=spin1+spin2+spin3+spin4
    U_delta_adim=2_sp*spin_value*val
end function U_delta_adim

! Magenitación
function M_adim(aux_matrix_pbc,n)
    use module_precision
    implicit none
    integer(sp), intent(in) :: n,aux_matrix_pbc(n+2,n+2)
    integer(sp)             :: i,j
    real(dp)                :: M_adim
    M_adim=0._dp
    do j=2,n;do i=2,n
        M_adim=M_adim+real(aux_matrix_pbc(i,j),dp)
    end do;end do
    ! M_adim=(1._dp/(n*n))*abs(M_adim) ! con valor absoluto
    M_adim=(1._dp/(n*n))*M_adim ! sin valor absoluto
end function M_adim

function M_exact_adim(n,T_adim,Tc_adim)
    use module_precision
    implicit none
    integer(sp), intent(in) :: n
    real(dp),    intent(in) :: T_adim,Tc_adim
    real(dp)                :: zz,M_exact_adim
    if (T_adim<Tc_adim) then
        zz=exp(-4*(1._dp/T_adim)) ! zz=z**2
        M_exact_adim=((1+zz)**0.25_dp)*((1-6*zz+zz*zz)**0.125_dp)*(1._dp/sqrt(1._dp-zz))
    else;M_exact_adim=0._dp;end if
end function M_exact_adim