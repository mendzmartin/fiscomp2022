!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
!
!
!----------------------------------------------------------

program EDO_primer_orden

	use module_presition 		! module for define presition
	use module_functions_1D 	! module for use 1D functions
	use module_functions_2D 	! module for use 2D functions
	use module_EDO_primer_orden	! module for resolve firste erder ODEs
	use module_numerical_error 	! module for estimate numerical errors

	implicit none
	
	! Explicit variables declaration
	integer(sp) 							:: n 				! points
	integer(sp) 							:: i,index_prev 		! index loop
	integer(sp)								:: function_type_EDO 	! tipo de función a usar en metodos numéricos
	integer(sp)								:: function_type_exact 	! tipo de función a usar en cálculo exacto
	integer(sp) 							:: istat				! integer of simple presition
	real(dp)								:: a,b					! intervalo de validez para la solución aproximada
	real(dp)								:: h					! step
	real(dp) 								:: y_0 					! initial condition
	real(dp) 								:: rel_error_euler 		! relative error Euler method
	real(dp) 								:: rel_error_RK2_hu		! relative error 2nd order RK-Hung method
	real(dp) 								:: rel_error_RK2_mp		! relative error 2nd order RK-Middle-Point method
	real(dp) 								:: rel_error_RK2_ra		! relative error 2nd order RK-Raltson method
	real(dp) 								:: rel_error_RK4_cl		! relative error 4th order RK-Classic method
	real(dp), dimension(:), allocatable 	:: y_exact 				! vector con solución exacta de la EDO
	real(dp), dimension(:), allocatable 	:: y_euler 				! aprox solution vector using Euler method
	real(dp), dimension(:), allocatable 	:: y_RK2_hu				! aprox solution vector using 2nd order RK-Hung method
	real(dp), dimension(:), allocatable 	:: y_RK2_mp				! aprox sol vector using 2nd order RK-Middle-Point method
	real(dp), dimension(:), allocatable 	:: y_RK2_ra				! aprox solution vector using 2nd order RK-Raltson method
	real(dp), dimension(:), allocatable 	:: y_RK4_cl 			! aprox soluction vector using 4th order RK-Classic method
	real(dp), dimension(:), allocatable 	:: x 					! variable vector
	
	! Datos pedidos al usuario
	write(*,*) 'Ingrese el limite inferior de integración (a) y presione Enter.'
	read(*,*) a
	write(*,*) 'Ingrese el limite superior de integración (a) y presione Enter.'
	read(*,*) b
	write(*,*) 'Ingrese la condición inicial dy/dx para x=0 y presione Enter.'
	read(*,*) y_0
	write(*,*) 'Ingrese la cantidad de puntos a utilizar (n) y presione Enter.'
	read(*,*) n

	
	
	open( 10, file = './result.dat', status = 'replace', iostat = istat )
	write(*,*) 'Input/Output file. istat = ', istat
	20 format (F10.2, F12.4, F12.4, F12.4, F12.4, F12.4, F12.4, E12.4, E12.4, E12.4, E12.4, E12.4)
	
	function_type_exact = 3 ! función 1D = exp(-x^2/2)
	function_type_EDO = 1 	! función 2D = -xy
	
	h = abs(b - a) * ( 1._dp / (n-1_sp) )

	allocate(x(n))
	x(1) = a
	allocate(y_exact(n))
	y_exact(1) = f_1D(x(1), function_type_exact)
	
	allocate(y_euler(n))
	call euler_method( n, a, b, y_0, y_euler, function_type_EDO )	! calculamos solución con método de euler
	allocate(y_RK2_hu(n))
	call RK2(1, n, a, b, y_0, y_RK2_hu, function_type_EDO) 			! calculamos solución con método de RK2
	allocate(y_RK2_mp(n))
	call RK2(2, n, a, b, y_0, y_RK2_mp, function_type_EDO) 			! calculamos solución con método de RK2
	allocate(y_RK2_ra(n))
	call RK2(3, n, a, b, y_0, y_RK2_ra, function_type_EDO) 			! calculamos solución con método de RK2
	allocate(y_RK4_cl(n))
	call RK4(1, n, a, b, y_0, y_RK4_cl, function_type_EDO)			! calculamos solución con método de RK4
	
	call basic_errors_num(y_exact(1), y_euler(1), rel_error_euler, 2) 		! calculamos error relativo método de euler
	call basic_errors_num(y_exact(1), y_RK2_hu(1), rel_error_RK2_hu, 2) 	! calculamos error relativo método de RK2
	call basic_errors_num(y_exact(1), y_RK2_mp(1), rel_error_RK2_mp, 2) 	! calculamos error relativo método de RK2
	call basic_errors_num(y_exact(1), y_RK2_ra(1), rel_error_RK2_ra, 2) 	! calculamos error relativo método de RK2
	call basic_errors_num(y_exact(1), y_RK4_cl(1), rel_error_RK4_cl, 2) 	! calculamos error relativo método de RK4
	
	write(10,20) 	x(1), y_exact(1),&
					y_euler(1),&
					y_RK2_hu(1), y_RK2_mp(1), y_RK2_ra(1),&
					y_RK4_cl(1),&
					rel_error_euler,&
					rel_error_RK2_hu, rel_error_RK2_mp, rel_error_RK2_ra,&
					rel_error_RK4_cl
	
	do i = 2, n, 1
		index_prev = i-1
		x(i) = x(index_prev) + h
		y_exact(i) = f_1D(x(i),function_type_exact)
		
		call basic_errors_num(y_exact(i), y_euler(i), rel_error_euler, 2) 		! calculamos error relativo método de euler
		call basic_errors_num(y_exact(i), y_RK2_hu(i), rel_error_RK2_hu, 2) 	! calculamos error relativo método de RK2
		call basic_errors_num(y_exact(i), y_RK2_mp(i), rel_error_RK2_mp, 2) 	! calculamos error relativo método de RK2
		call basic_errors_num(y_exact(i), y_RK2_ra(i), rel_error_RK2_ra, 2) 	! calculamos error relativo método de RK2
		call basic_errors_num(y_exact(i), y_RK4_cl(i), rel_error_RK4_cl, 2) 	! calculamos error relativo método de RK4
		
		write(10,20) 	x(i), y_exact(i),&
						y_euler(i),&
						y_RK2_hu(i), y_RK2_mp(i), y_RK2_ra(i),&
						y_RK4_cl(i),&
						rel_error_euler,&
						rel_error_RK2_hu, rel_error_RK2_mp, rel_error_RK2_ra,&
						rel_error_RK4_cl
	end do
	
	deallocate(y_exact, x)
	deallocate(y_euler)
	deallocate(y_RK2_hu, y_RK2_mp, y_RK2_ra)
	deallocate(y_RK4_cl)
	
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
! gfortran -g -fcheck=all -Wall -o EDO_primer_orden ../../modules/module_presition.f90 ../../modules/module_functions_1D.f90 ../../modules/module_functions_2D.f90 ../../modules/module_EDO_primer_orden.f90 ../../modules/module_numerical_error.f90 EDO_primer_orden.f90 && ./EDO_primer_orden

! gfortran -o EDO_primer_orden ../../modules/module_presition.f90 ../../modules/module_functions_1D.f90 ../../modules/module_functions_2D.f90 ../../modules/module_EDO_primer_orden.f90 ../../modules/module_numerical_error.f90 EDO_primer_orden.f90 && ./EDO_primer_orden
!----------------------------------------------------------
