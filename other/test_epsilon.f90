!----------------------------------------------------------
! epsilon(x) returns the smallest number E of the same kind
!  as x such that 1 + E > 1.
!----------------------------------------------------------
program test_epsilon

    use module_presition

    implicit none
    
    real(sp) :: value_sp = 2.357 ! single presition
    real(dp) :: value_dp = 2.357 ! double presition
    
    20 format (A20, E11.3)
    
    write(*,20) "epsilon_machine_sp = ", epsilon(value_sp)
    write(*,20) "epsilon_machine_dp = ", epsilon(value_dp)
end program test_epsilon

!----------------------------------------------------------
! Notas
!----------------------------------------------------------
! link = https://gcc.gnu.org/onlinedocs/gfortran/EPSILON.html

!----------------------------------------------------------
! Regla de compilación
!----------------------------------------------------------
! gfortran -o test_epsilon.o test_epsilon.f90 ../modules/module_presition.f90 && ./test_epsilon.o
