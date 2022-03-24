!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
!
!
!----------------------------------------------------------

program EDO_primer_orden

	use module_presition 		! use module of presition
	use module_functions 		! use module of functions

	implicit none
	
	! Explicit variables declaration
	real(dp)				:: a,b					! intervalo de validez para la solución aproximada
	real(dp)				:: h					! step
	integer(sp), parameter 	:: n = 10 				! points
	real(dp), dimension(n) 	:: y_exact 				! vector con solución exacta de la EDO
	real(dp), dimension(n) 	:: y_euler 				! vector con solución aproximada de la EDO usando método de Euler
	real(dp), dimension(n) 	:: x 					! variable vector
	real(dp) 				:: y_0 					! initial condition
	integer(sp) 			:: i,index_prev 		! index loop
	integer(sp)				:: function_type_euler 	! tipo de función a usar en método de Euler
	integer(sp)				:: function_type_exact 	! tipo de función a usar en cálculo exacto
	integer(sp) 			:: istat				! integer of simple presition
	
	! Datos pedidos al usuario
	write( *, * ) 'Ingrese el limite inferior de integración (a) y presione Enter.'
	read( *, * ) a
	write( *, * ) 'Ingrese el limite superior de integración (a) y presione Enter.'
	read( *, * ) b
	write( *, * ) 'Ingrese la condición inicial dy/dx para x=0 y presione Enter.'
	read( *, * ) y_0
	
	function_type_exact = 3 ! exp(-x^2/2)
	function_type_euler = 1 ! -xy
	
	h = ( b - a ) * ( 1._dp / n )
	
	x(1) = a
	y_euler(1) = y_0
	y_exact(1) = f(x(1), function_type_exact)
	
	open( 10, file = './result.dat', status = 'replace', iostat = istat )
	write(*,*) 'Input/Output file. istat = ', istat
	20 format (E15.6, E15.6, E15.6)
	
	write(10,20) x(1), y_exact(1), y_euler(1)
	
	do i = 2, n, 1
		
		index_prev = i-1
		
		x(i) = x(index_prev) + h
		
		y_exact(i) = f(x(i),function_type_exact)
		
		y_euler(i) = y_euler(index_prev) + h*g(x(index_prev), y_euler(index_prev), function_type_euler)
		
		write(10,20) x(i), y_exact(i), y_euler(i)
		
	end do
	
	close(10)
	
end program EDO_primer_orden



!----------------------------------------------------------
! REFERENCES
!----------------------------------------------------------
!
!
!
!----------------------------------------------------------

!----------------------------------------------------------
! COMPILATION RUL
!----------------------------------------------------------
! gfortran -o EDO_primer_orden ../../modules/module_presition.f90 ../../modules/module_functions.f90 EDO_primer_orden.f90 && ./EDO_primer_orden
!
!
!----------------------------------------------------------
