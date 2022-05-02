program test_module_random_generator
    use module_precision
    use module_random_generator
    use module_mzran
    use module_mt19937
    implicit none
    integer(sp) :: seed
    seed=121648
    ! test ran0
    19 format(A10,I10)
    20 format(A10,E15.8)
    write(*,19) 'seed=',seed
    write(*,20) 'ran0=',ran0(seed)
    write(*,20) 'ran2=',ran2(seed)
    write(*,20) 'mzran=',rmzran()
    call sgrnd(seed); write(*,20) 'mt=',grnd()
    
end program test_module_random_generator
! RUN
! gfortran -o test_module_random_generator.o ../module_precision.f90 ../module_random_generator.f90 ../module_mzran.f90 ../module_mt19937.f90 test_module_random_generator.f90
! ./test_module_random_generator.o && rm *.mod *.o