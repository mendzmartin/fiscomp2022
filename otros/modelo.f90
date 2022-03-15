!----------------------------------------------------------
! Programa modelo para entender como se deben declarar
!  los parámetros y como se usan variables lógicas
!  booleanas.
!----------------------------------------------------------

program modelo
	implicit none
	integer,parameter	:: np = kind( 1.0 ) ! np=4 es Int32 precisión simple
	logical				:: p1, p2
	real(np)			:: a, b, c
	integer				:: n, m
	p1 = n**m < a*b
	p1 = ( ( a**n ) == ( b**m ) )
	write(*,*) p1, p2
end program modelo

!----------------------------------------------------------
! Regla de compilación
!----------------------------------------------------------
! gfortran -o modelo modelo.f90 && ./modelo
