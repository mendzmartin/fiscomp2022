!----------------------------------------------------------
! epsilon(x) returns the smallest number E of the same kind
!  as x such that 1 + E > 1.
!----------------------------------------------------------
program test_epsilon
    real	:: x = 3.143
    real(8) :: y = 2.33
    print *, epsilon(x)
    print *, epsilon(y)
end program test_epsilon

!----------------------------------------------------------
! Notas
!----------------------------------------------------------
! link = https://gcc.gnu.org/onlinedocs/gfortran/EPSILON.html

!----------------------------------------------------------
! Regla de compilación
!----------------------------------------------------------
! gfortran -o modelo modelo.f90 && ./modelo
