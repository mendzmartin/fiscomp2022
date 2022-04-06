!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
!
!
!----------------------------------------------------------

program EDO_primer_orden_error_analysis

	use module_presition 		! module for define presition
	use module_functions_1D 	! module for use 1D functions
	use module_functions_2D 	! module for use 2D functions
	use module_EDO_primer_orden	! module for resolve firste erder ODEs
	use module_numerical_error 	! module for estimate numerical errors

	implicit none
	
	! Explicit variables declaration
	
	! number of points to obtain an absolute error [%] in the order of (10^-5)
	integer(sp), parameter 	::  n_euler = 35000 ! used points in Euler method
	integer(sp), parameter 	:: n_RK2_hu = 200 	! used points in 2nd order RK-Hung method
	integer(sp), parameter 	:: n_RK2_mp = 250 	! used points in 2nd order RK-Middle-Point method
	integer(sp), parameter 	:: n_RK2_ra = 50 	! used points in 2nd order RK-Raltson method
	integer(sp), parameter 	:: n_RK4_cl = 4 	! used points in 4th order RK-Classic method

	real(dp) :: a, b			! intervalo de validez para la solución aproximada
	real(dp) :: y_0 			! initial condition

	
	! Datos pedidos al usuario
	write( *, * ) 'Ingrese el limite inferior de integración (a) y presione Enter.'
	read( *, * ) a
	write( *, * ) 'Ingrese el limite superior de integración (a) y presione Enter.'
	read( *, * ) b
	write( *, * ) 'Ingrese la condición inicial dy/dx para x=0 y presione Enter.'
	read( *, * ) y_0

	! Analysis of Euler method
	call results_collection(10, './result_euler.dat', n_euler, a, b, y_0, 1)
	
	! Analysis of 2nd order RK-Hung method
	call results_collection(11, './result_RK2_hu.dat', n_RK2_hu, a, b, y_0, 2)

	! Analysis of 2nd order RK-Middle-Point method
	call results_collection(12, './result_RK2_mp.dat', n_RK2_mp, a, b, y_0, 3)

	! Analysis of 2nd order RK-Raltson method
	call results_collection(13, './result_RK2_ra.dat', n_RK2_ra, a, b, y_0, 4)

	! Analysis of 4th order RK-Classic method
	call results_collection(14, './result_RK4_cl.dat', n_RK4_cl, a, b, y_0, 5)

end program EDO_primer_orden_error_analysis


subroutine results_collection(unit_number, file_name, n, a, b, y_0, method_type)

	use module_presition 		! module for define presition
	use module_functions_1D 	! module for use 1D functions
	use module_functions_2D 	! module for use 2D functions
	use module_EDO_primer_orden	! module for resolve firste erder ODEs
	use module_numerical_error 	! module for estimate numerical errors
	
	implicit none

	! Data dictionary: declare calling parameter types & definitions
	integer(sp), 	intent(in) :: unit_number
	character(20), 	intent(in) :: file_name
	integer(sp), 	intent(in) :: n 			! used points in numercal method
	integer(sp), 	intent(in) :: method_type 	! type of used method
	real(dp), 		intent(in) :: a,b 			! intervalo de validez para la solución aproximada
	real(dp), 		intent(in) :: y_0 			! initial condition

	! Data dictionary: declare local variables types & definitions
	integer(sp) 			:: i,index_prev 			! index loop
	integer(sp)				:: function_type_EDO 		! tipo de función a usar en metodos numéricos
	integer(sp)				:: function_type_exact 		! tipo de función a usar en cálculo exacto
	integer(sp) 			:: istat					! integer of simple presition
	real(dp)				:: h						! step
	real(dp) 				:: rel_error 				! relative error
	real(dp) 				:: start_time, finish_time 	! variables to compute the elapsed CPU time in seconds
	real(dp) 				:: elapsed_time 			! elapsed CPU time
	real(dp), dimension(n) 	:: y_method 				! vector con solución aproximada de la EDO
	real(dp), dimension(n) 	:: y_exact 					! vector con solución exacta de la EDO
	real(dp), dimension(n) 	:: x						! variable vector
	
	
	function_type_exact = 3 ! función 1D = exp(-x^2/2)
	function_type_EDO = 1 	! función 2D = -xy
	
	open( unit = unit_number, file = file_name, status = 'replace', iostat = istat )
	write(*,*) 'Input/Output ', unit_number, 'file. istat = ', istat
	20 format (F10.2, F12.4, F12.4, E12.4)
	21 format (F10.2, F12.4, F12.4, E12.4, E12.4)
	
	x(1) = a
	y_exact(1) = f_1D(x(1), function_type_exact)

	select case(method_type)
		case(1)
			call cpu_time(start_time)
			call euler_method( n, a, b, y_0, y_method, function_type_EDO )
			call cpu_time(finish_time)
		case(2)
			call cpu_time(start_time)
			call RK2(1, n, a, b, y_0, y_method, function_type_EDO)
			call cpu_time(finish_time)
		case(3)
			call cpu_time(start_time)
			call RK2(2, n, a, b, y_0, y_method, function_type_EDO)
			call cpu_time(finish_time)
		case(4)
			call cpu_time(start_time)
			call RK2(3, n, a, b, y_0, y_method, function_type_EDO)
			call cpu_time(finish_time)
		case(5)
			call cpu_time(start_time)
			call RK4(1, n, a, b, y_0, y_method, function_type_EDO)
			call cpu_time(finish_time)
		case default
			write(*,*) 'Invalid method type'
		end select
	
	elapsed_time = 0.001_dp*(finish_time - start_time) ! elapsed CPU time in mili-seconds
	
	call basic_errors_num(y_exact(1), y_method(1), rel_error, 2)
	
	write(unit_number,21) 	x(1), y_exact(1), y_method(1), rel_error, elapsed_time
	
	h = abs(b - a) * ( 1._dp / (n-1_sp) )
	
	do i = 2, n, 1
		index_prev = i-1
		x(i) = x(index_prev) + h
		y_exact(i) = f_1D(x(i),function_type_exact)
		
		call basic_errors_num(y_exact(i), y_method(i), rel_error, 2)
		
		write(unit_number,20) 	x(i), y_exact(i), y_method(i), rel_error
	end do
	
	close(unit_number)

end subroutine results_collection



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
! gfortran -g -fcheck=all -Wall -o EDO_primer_orden_error_analysis ../../modules/module_presition.f90 ../../modules/module_functions_1D.f90 ../../modules/module_functions_2D.f90 ../../modules/module_EDO_primer_orden.f90 ../../modules/module_numerical_error.f90 EDO_primer_orden_error_analysis.f90 && ./EDO_primer_orden_error_analysis

! gfortran -o EDO_primer_orden_error_analysis ../../modules/module_presition.f90 ../../modules/module_functions_1D.f90 ../../modules/module_functions_2D.f90 ../../modules/module_EDO_primer_orden.f90 ../../modules/module_numerical_error.f90 EDO_primer_orden_error_analysis.f90 && ./EDO_primer_orden_error_analysis
!----------------------------------------------------------
