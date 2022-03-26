module module_EDO_primer_orden

	use module_presition
	use module_functions_2D
	
	implicit none
	
	contains
	
	subroutine euler_method ( n, a, b, y_0, y_euler, function_type )
		
		! Data dictionary: declare calling parameter types & definitions
		integer(sp), intent(in)					:: n 				! points
		integer(sp), intent(in) 				:: function_type 	! type of function to resolve
		real(dp), intent(in)	 				:: y_0 				! initial condition
		real(dp), intent(in)					:: a,b				! validity range of approximate solution
		real(dp), dimension(n), intent(out) 	:: y_euler 			! aproximate solution vector of ODE using Euler method
		
		! Data dictionary: declare local variables types & definitions
		integer(sp) 							:: i,index_prev 	! index loop
		real(dp)								:: h				! step
		real(dp), dimension(n) 					:: x 				! grid points vector
		
		! Execution section
		
		h = abs(b - a) * ( 1._dp / (n-1_sp) )
		
		x(1) = a
		y_euler(1) = y_0
		
		do i = 2, n, 1
			index_prev = i-1
			x(i) = x(index_prev) + h
			
			! ¡¡OJO!! se usa la f(x,y) porque se sabe de antemano que hay que usar esta
			!  function del modulo "module_functions_2D" si se quiere trabajar con una
			!  f(x) se debe agregar el modulo "module_functions_1D"
			
			y_euler(i) = y_euler(index_prev) + h*f_2D(x(index_prev), y_euler(index_prev), function_type)
		end do
		
	end subroutine euler_method
	
	subroutine RK2(RK2_type, n, a, b, y_0, y_RK2, function_type)

		! Data dictionary: declare calling parameter types & definitions
		integer(sp), 			intent(in)		:: RK2_type
		integer(sp), 			intent(in)		:: n 				! points
		integer(sp), 			intent(in) 		:: function_type 	! tipo de función a integrar
		real(dp), 				intent(in)		:: y_0 				! initial condition
		real(dp), 				intent(in)		:: a,b				! intervalo de validez para la solución aproximada
		real(dp), dimension(n), intent(out) 	:: y_RK2 			! vector con solución aproximada de la EDO usando el método de Ruge Kutta de 2do orden 
		
		! Data dictionary: declare local variables types & definitions
		integer(sp) 			:: i,index_prev 	! index loop
		real(dp) 				:: k_1, k_2 		! recurrence relations
		real(dp) 				:: y_improved 		! improved function evaluation
		real(dp) 				:: h_improved 		! improved step evaluation
		real(dp), dimension(n) 	:: x 				! grid points vector
		real(dp) 				:: h 				! step
		
	
		select case (RK2_type)
			case(1) ! Método de Heun con un solo corrector (a_2 = 1/2)
				h = abs(b - a) * ( 1._dp / (n-1_sp) )
				x(1) = a
				y_RK2(1) = y_0
			
				do i = 2, n, 1
					index_prev = i-1
					x(i) = x(index_prev) + h
					
					! ¡¡OJO!! se usa la f(x,y) porque se sabe de antemano que hay que usar esta
					!  function del modulo "module_functions_2D" si se quiere trabajar con una
					!  f(x) se debe agregar el modulo "module_functions_1D"
					
					k_1 = f_2D(x(index_prev), y_RK2(index_prev), function_type)
					y_improved = y_RK2(index_prev) + (k_1*h)
					k_2 = f_2D(x(i), y_improved, function_type )
					
					! y(i+1) = (y(i) + increment_function)
					y_RK2(i) = y_RK2(index_prev) + (0.5_dp)*h*(k_1 + k_2)
				end do
				
			case(2) ! Método del punto medio (a_2 = 1)
				h = abs(b - a) * ( 1._dp / (n-1_sp) )
				x(1) = a
				y_RK2(1) = y_0
			
				do i = 2, n, 1
					index_prev = i-1
					x(i) = x(index_prev) + h
					
					! ¡¡OJO!! se usa la f(x,y) porque se sabe de antemano que hay que usar esta
					!  function del modulo "module_functions_2D" si se quiere trabajar con una
					!  f(x) se debe agregar el modulo "module_functions_1D"
					
					k_1 = f_2D(x(index_prev), y_RK2(index_prev), function_type)
					h_improved = 0.5_dp*h
					y_improved = y_RK2(index_prev) + (k_1*h_improved)
					k_2 = f_2D( (x(i)-h_improved), y_improved, function_type )
					
					! y(i+1) = (y(i) + increment_function)
					y_RK2(i) = y_RK2(index_prev) + h*k_2
				end do
				
			case(3) ! Método de Ralston (a_2 = 2/3)
				h = abs(b - a) * ( 1._dp / (n-1_sp) )
				x(1) = a
				y_RK2(1) = y_0
			
				do i = 2, n, 1
					index_prev = i-1
					x(i) = x(index_prev) + h
					
					! ¡¡OJO!! se usa la f(x,y) porque se sabe de antemano que hay que usar esta
					!  function del modulo "module_functions_2D" si se quiere trabajar con una
					!  f(x) se debe agregar el modulo "module_functions_1D"
					
					k_1 = f_2D(x(index_prev), y_RK2(index_prev), function_type)
					h_improved = 0.25_dp*h
					y_improved = y_RK2(index_prev) + (3._dp*k_1*h_improved)
					k_2 = f_2D( (x(i)-h_improved), y_improved, function_type )
					
					! y(i+1) = (y(i) + increment_function)
					y_RK2(i) = y_RK2(index_prev) + (1._dp/3._dp)*h*(k_1+2*k_2)
				end do
				
			case default
				write(*,*) 'Invalid RK2 type'
		end select
	end subroutine RK2
	
	subroutine RK4(RK4_type, n, a, b, y_0, y_RK4, function_type)
		! Data dictionary: declare calling parameter types & definitions
		integer(sp), 			intent(in)		:: RK4_type
		integer(sp), 			intent(in)		:: n 				! points
		integer(sp), 			intent(in) 		:: function_type 	! tipo de función a integrar
		real(dp), 				intent(in)		:: y_0 				! initial condition
		real(dp), 				intent(in)		:: a,b				! intervalo de validez para la solución aproximada
		real(dp), dimension(n), intent(out) 	:: y_RK4 			! vector con solución aproximada de la EDO usando el método de Ruge Kutta de 2do orden 
		
		! Data dictionary: declare local variables types & definitions
		integer(sp) 			:: i,index_prev 		! index loop
		real(dp) 				:: k_1, k_2, k_3, k_4	! recurrence relations
		real(dp) 				:: y_improved_k2		! improved function evaluation
		real(dp) 				:: y_improved_k3		! improved function evaluation
		real(dp) 				:: y_improved_k4		! improved function evaluation
		real(dp) 				:: h_improved 			! improved step evaluation
		real(dp), dimension(n) 	:: x 					! grid points vector
		real(dp) 				:: h 					! step
		
		select case (RK4_type)
			case(1) ! Método clásico
				h = abs(b - a) * ( 1._dp / (n-1_sp) )
				x(1) = a
				y_RK4(1) = y_0
			
				do i = 2, n, 1
					index_prev = i-1
					x(i) = x(index_prev) + h
					
					! ¡¡OJO!! se usa la f(x,y) porque se sabe de antemano que hay que usar esta
					!  function del modulo "module_functions_2D" si se quiere trabajar con una f(x)
					!  se debe agregar el modulo "module_functions_1D" y cambiar f_2D por f_1D
					
					k_1 = f_2D(x(index_prev), y_RK4(index_prev), function_type)
					
					h_improved = 0.5_dp*h
					
					y_improved_k2 = y_RK4(index_prev) + (k_1*h_improved)
					k_2 = f_2D(x(i)-h_improved, y_improved_k2, function_type )
					
					y_improved_k3 = y_RK4(index_prev) + (k_2*h_improved)
					k_3 = f_2D(x(i)-h_improved, y_improved_k3, function_type )
					
					y_improved_k4 = y_RK4(index_prev) + (k_3*h)
					k_4 = f_2D(x(i), y_improved_k4, function_type )
					
					! y(i+1) = (y(i) + increment_function)
					y_RK4(i) = y_RK4(index_prev) + (1._dp/6._dp)*h*(k_1 + 2._dp*(k_2+k_3)+k_4)
				end do
				
			case default
				write(*,*) 'Invalid RK4 type'
		
		end select
				
	end subroutine RK4
	
end module module_EDO_primer_orden
