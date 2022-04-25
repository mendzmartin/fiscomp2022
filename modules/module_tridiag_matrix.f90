module module_tridiag_matrix
    use module_precision
    implicit none
    contains
    subroutine implicit_method(n,diag,diag_sup,diag_inf,data_vector,unknown_vector)
        ! variables de entrada/salida
        integer(sp),    intent(in)  :: n                               ! dimension del vector diagonal
        real(dp),       intent(out) :: unknown_vector(n)               ! vector de incognitas
        real(dp),       intent(in)  :: diag(n),diag_sup(n),diag_inf(n) ! diagonales central, superior e inferior
        real(dp),       intent(in)  :: data_vector(n)                  ! vector de datos conocidos
        ! variables locales
        integer(sp) :: i                    ! variable del loop
        real(dp)    :: diag_sup_new(n)      ! diagonal superior redefinida
        real(dp)    :: data_vector_new(n)   ! vector de datos redefinidos
        real(dp)    :: factor
        ! descomposición LU
        factor=1._dp/diag(1)
        diag_sup_new(1)=diag_sup(1)*factor
        data_vector_new(1)=data_vector(1)*factor
        do i=2,n
            factor=1._dp/(diag(i)-diag_inf(i)*diag_sup_new(i-1))
            diag_sup_new(i)=diag_sup(i)*factor
            data_vector_new(i)=(data_vector(i)-diag_inf(i)*data_vector_new(i-1))*factor
        end do
        ! sustitución hacia atrás
        unknown_vector(n)=data_vector_new(n)
        do i=(n-1),1,-1
            unknown_vector(i)=data_vector_new(i)-diag_sup_new(i)*unknown_vector(i+1)
        end do
    end subroutine implicit_method
end module module_tridiag_matrix
! Los vectores de entrada deben definirse de la siguiente manera
! diag_sup = [ ds(1)   ds(2) ... ds(n-1) ds(n)=0 ]
! diag_inf = [ di(1)=0 di(2) ... di(n-1) di(n)   ]
! diag     = [ d(1)    d(2)  ... di(n-1) d(n)    ]