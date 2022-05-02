program test_module_tridiag_matrix
    use module_precision
    use module_tridiag_matrix
    implicit none
    real(dp) :: diag(3),diag_inf(3),diag_sup(3)
    real(dp) :: vector1(3),vector2(3)
    diag=[2,3,5]
    diag_sup=[7,11,0]
    diag_inf=[0,17,13]
    vector1=[1,2,3]
    call implicit_method(3,diag,diag_sup,diag_inf,vector1,vector2)
    write(*,*) vector1
    write(*,*) vector2
end program test_module_tridiag_matrix
! run command
! gfortran -o test_module_tridiag_matrix.o module_precision.f90 module_tridiag_matrix.f90 test_module_tridiag_matrix.f90 && ./test_module_tridiag_matrix.o && rm *.mod *.o