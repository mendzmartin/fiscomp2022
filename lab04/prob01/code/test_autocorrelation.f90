! Test para la funcion de Autocorrelación
! make clean && make test_autocorrelation.o && ./test_autocorrelation.o
program test_autocorrelation
    use module_precision
    implicit none

    integer(sp), parameter :: num_of_terms=3_sp,tau_corr=3_sp,m_exp=2_sp
    real(dp)               :: obs
    integer(sp)            :: i,j
    real(dp)               :: autocor_vector_mom(tau_corr)
    real(dp)               :: autocor_vector(tau_corr)
    real(dp)               :: aux_vector1(tau_corr),aux_vector2(tau_corr)
    real(dp)               :: mask_vector(tau_corr)

    autocor_vector_mom(:)=0._dp

    do j=1,m_exp ! cantidad de experimentos
        do i=1,(num_of_terms+tau_corr)
            obs=real(i,dp)
            call func_autocor(obs,i,tau_corr,num_of_terms,autocor_vector,aux_vector1,mask_vector,aux_vector2)
        end do

        ! do i=1,tau_corr
        !     write(*,*) i,autocor_vector(i)
        ! end do
        ! el resultado debería ser autocor_vector=[70/5,50/4,32/3]

        autocor_vector_mom(:)=autocor_vector_mom(:)+autocor_vector(:)
    end do
    autocor_vector_mom(:)=autocor_vector_mom(:)*(1._dp/real(m_exp,dp))

    do i=1,tau_corr
        write(*,*) i,autocor_vector_mom(i)
    end do

end program test_autocorrelation

subroutine func_autocor(obs,time_index,tau_corr,num_of_terms,&
                        autocor_vector,aux_vector1,mask_vector,aux_vector2)
    use module_precision
    implicit none

    integer(sp), intent(in)    :: time_index ! indice del observable
    integer(sp), intent(in)    :: num_of_terms
    integer(sp), intent(in)    :: tau_corr   ! tiempo de correlación
    real(dp),    intent(in)    :: obs        ! observable
    real(dp),    intent(inout) :: autocor_vector(tau_corr)
    real(dp),    intent(inout) :: aux_vector1(tau_corr),mask_vector(tau_corr)
    real(dp),    intent(inout) :: aux_vector2(tau_corr)
    integer(sp)   :: i,index

    ! llenamos por primera vez el vector de dim=tau_corr
    if (time_index<=tau_corr) autocor_vector(time_index)=obs
    if (time_index==tau_corr) then 
        aux_vector1(:)=autocor_vector(:)
        aux_vector2(:)=autocor_vector(:)
    end if

    ! computamos las correlaciones a todo tiempo
    if (time_index>tau_corr.and.time_index<=(num_of_terms+tau_corr)) then
        do i=1_sp,tau_corr-1_sp;if (mod(time_index,tau_corr)==i) index=i;end do
        if (mod(time_index,tau_corr)==0_sp) index=tau_corr
        aux_vector2(:)=cshift(aux_vector2(:),shift=1)
        aux_vector2(tau_corr)=obs
        autocor_vector(index)=dot_product(aux_vector1(:),aux_vector2(:))
    end if

    ! computamos las correlaciones faltantes
    if (time_index==(num_of_terms+tau_corr)) then
        mask_vector(:)=1._dp
        aux_vector1(:)=aux_vector2(:)
        do index=1,tau_corr-1
            aux_vector2(:)=cshift(aux_vector2(:),shift=1)
            mask_vector(tau_corr-(index-1))=0._dp
            aux_vector2(:)=aux_vector2(:)*mask_vector(:)
            autocor_vector(index)=autocor_vector(index)+dot_product(aux_vector1(:),aux_vector2(:))
            autocor_vector(index)=autocor_vector(index)*(1._dp/real(num_of_terms+tau_corr-index,dp))
        end do
        autocor_vector(tau_corr)=autocor_vector(tau_corr)*(1._dp/real(num_of_terms,dp))
    end if

end subroutine func_autocor