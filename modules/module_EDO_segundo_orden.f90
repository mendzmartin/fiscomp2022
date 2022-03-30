module module_EDO_segundo_orden

	use module_presition
	use module_functions_2D
	
	implicit none
	
	contains
	
	subroutine euler_method( n, a, b, y1_0, y2_0, y1_euler, y2_euler, function_type, input_type )
		
		! Data dictionary: declare calling parameter types & definitions
		integer(sp), intent(in)					:: n 					! points
		integer(sp), intent(in) 				:: function_type 		! type of function to resolve
		integer(sp), intent(in) 				:: input_type 			! type of input to create function
		real(dp), intent(in)	 				:: y1_0, y2_0 			! initial condition
		real(dp), intent(in)					:: a,b					! validity range of approximate solution
		real(dp), dimension(n), intent(out) 	:: y1_euler, y2_euler	! aproximate solution vector of ODE using Euler method
		
		
		! Data dictionary: declare local variables types & definitions
		integer(sp) 							:: i,index_prev 	! index loop
		real(dp)								:: h				! step
		real(dp), dimension(n) 					:: x 				! grid points vector
		
		! Execution section
		
		h = abs(b - a) * ( 1._dp / (n-1_sp) )
		
		x(1) = a
		y1_euler(1) = y1_0
		y2_euler(1) = y2_0
		
	select case(input_type)
		case(1)
			do i = 2, n, 1
				index_prev = i-1
				x(i) = x(index_prev) + h
				y1_euler(i) = y1_euler(index_prev) + h*(-y2_euler(index_prev))
				y2_euler(i) = y2_euler(index_prev) + h*y1_euler(index_prev)
			end do
		!case(2) ! using module_functions_2D
		!case(3) ! using module_functions_1D
		case default
			write(*,*) 'Invalid input type'
	end select
		
	end subroutine euler_method
	
	subroutine RK2(RK2_type, n, a, b, y1_0, y2_0, y1_RK2, y2_RK2, function_type, input_type)

		! Data dictionary: declare calling parameter types & definitions
		integer(sp), 			intent(in)		:: RK2_type
		integer(sp), 			intent(in)		:: n 				! points
		integer(sp), 			intent(in) 		:: function_type 	! tipo de función a integrar
		integer(sp), 			intent(in)		:: input_type 		! type of input to create function
		real(dp), 				intent(in)		:: y1_0 			! initial condition
		real(dp), 				intent(in)		:: y2_0 			! initial condition
		real(dp), 				intent(in)		:: a,b				! intervalo de validez para la solución aproximada
		real(dp), dimension(n), intent(out) 	:: y1_RK2 			! vector con solución aproximada
		real(dp), dimension(n), intent(out) 	:: y2_RK2 			! vector con solución aproximada 
		
		! Data dictionary: declare local variables types & definitions
		integer(sp) 			:: i,index_prev 	! index loop
		real(dp) 				:: k1_1, k1_2 		! recurrence relations
		real(dp) 				:: k2_1, k2_2 		! recurrence relations
		real(dp) 				:: y1_improved 		! improved function evaluation
		real(dp) 				:: y2_improved 		! improved function evaluation
		real(dp) 				:: h_improved 		! improved step evaluation
		real(dp), dimension(n) 	:: x 				! grid points vector
		real(dp) 				:: h 				! step
		
		select case(input_type)
			case(1) ! custom functions
				select case (RK2_type)
					case(1) ! Método de Heun con un solo corrector (a_2 = 1/2)
						h = abs(b - a) * ( 1._dp / (n-1_sp) )
						x(1) = a
						y1_RK2(1) = y1_0
						y2_RK2(1) = y1_0
					
						do i = 2, n, 1
							index_prev = i-1
							x(i) = x(index_prev) + h
							
							! ki_1 = fi(y)
							k1_1 = -y2_RK2(index_prev) 	! custom function f1(y)
							k2_1 = y1_RK2(index_prev)	! custon function f2(y)
							
							! yi_improved = yi_RK2(index_prev) + ki_1*h
							y1_improved = y1_RK2(index_prev) + (k1_1*h)
							y2_improved = y2_RK2(index_prev) + (k2_1*h)
							
							! ki_2 = fi(y_improve)
							k1_2 = -y2_improved
							k2_2 = y1_improved
							
							! y(i+1) = (y(i) + increment_function)
							y1_RK2(i) = y1_RK2(index_prev) + (0.5_dp)*h*(k2_1 + k2_2)
							y2_RK2(i) = y2_RK2(index_prev) + (0.5_dp)*h*(k2_1 + k2_2)
						end do
					case(2) ! Método del punto medio (a_2 = 1)
						h = abs(b - a) * ( 1._dp / (n-1_sp) )
						x(1) = a
						y1_RK2(1) = y1_0
						y2_RK2(1) = y2_0
					
						do i = 2, n, 1
							index_prev = i-1
							x(i) = x(index_prev) + h
							
							k1_1 = -y2_RK2(index_prev)
							k2_1 = y1_RK2(index_prev)
							h_improved = 0.5_dp*h
							y1_improved = y1_RK2(index_prev) + (k1_1*h_improved)
							y2_improved = y2_RK2(index_prev) + (k2_1*h_improved)
							k1_2 = -y2_improved
							k2_2 = y1_improved
							
							! y(i+1) = (y(i) + increment_function)
							y1_RK2(i) = y1_RK2(index_prev) + h*k1_2
							y2_RK2(i) = y2_RK2(index_prev) + h*k2_2
						end do
					case(3) ! Método de Ralston (a_2 = 2/3)
						h = abs(b - a) * ( 1._dp / (n-1_sp) )
						x(1) = a
						y1_RK2(1) = y1_0
						y2_RK2(1) = y2_0
					
						do i = 2, n, 1
							index_prev = i-1
							x(i) = x(index_prev) + h
							
							k1_1 = -y2_RK2(index_prev)
							k2_1 = y1_RK2(index_prev)
							h_improved = 0.25_dp*h
							y1_improved = y1_RK2(index_prev) + (3._dp*k1_1*h_improved)
							y2_improved = y2_RK2(index_prev) + (3._dp*k2_1*h_improved)
							k1_2 = -y2_improved
							k2_1 = y1_improved
							
							! y(i+1) = (y(i) + increment_function)
							y1_RK2(i) = y1_RK2(index_prev) + (1._dp/3._dp)*h*(k1_1+2*k1_2)
							y2_RK2(i) = y2_RK2(index_prev) + (1._dp/3._dp)*h*(k2_1+2*k2_2)
						end do
					case default
						write(*,*) 'Invalid RK2 type'
				end select
			! agregar casos si se quiere usar module_functions_1D o module_functions_2D
			!case(2) 
			case default
				write(*,*) 'Invalid input type'
		end select
	end subroutine RK2
	
	subroutine RK4(RK4_type, n, a, b, y1_0, y2_0, y1_RK4, y2_RK4, function_type, input_type)
		! Data dictionary: declare calling parameter types & definitions
		integer(sp), 			intent(in)		:: RK4_type
		integer(sp), 			intent(in)		:: n 				! points
		integer(sp), 			intent(in) 		:: function_type 	! tipo de función a integrar
		integer(sp), 			intent(in)		:: input_type 		! type of input to create function
		real(dp), 				intent(in)		:: y1_0 			! initial condition
		real(dp), 				intent(in)		:: y2_0 			! initial condition
		real(dp), 				intent(in)		:: a,b				! intervalo de validez para la solución aproximada
		real(dp), dimension(n), intent(out) 	:: y1_RK4 			! vector con solución aproximada
		real(dp), dimension(n), intent(out) 	:: y2_RK4 			! vector con solución aproximada
		
		! Data dictionary: declare local variables types & definitions
		integer(sp) 			:: i,index_prev 			! index loop
		real(dp) 				:: k1_1, k1_2, k1_3, k1_4	! recurrence relations
		real(dp) 				:: k2_1, k2_2, k2_3, k2_4	! recurrence relations
		real(dp) 				:: y1_improved_k2			! improved function evaluation
		real(dp) 				:: y1_improved_k3			! improved function evaluation
		real(dp) 				:: y1_improved_k4			! improved function evaluation
		real(dp) 				:: y2_improved_k2			! improved function evaluation
		real(dp) 				:: y2_improved_k3			! improved function evaluation
		real(dp) 				:: y2_improved_k4			! improved function evaluation
		real(dp) 				:: h_improved 				! improved step evaluation
		real(dp), dimension(n) 	:: x 						! grid points vector
		real(dp) 				:: h 						! step
		
		select case(input_type)
			case(1) ! custom functions
				select case (RK4_type)
					case(1) ! Método clásico
						h = abs(b - a) * ( 1._dp / (n-1_sp) )
						x(1) = a
						y1_RK4(1) = y1_0
						y2_RK4(1) = y2_0
					
						do i = 2, n, 1
							index_prev = i-1
							x(i) = x(index_prev) + h
							
							k1_1 = -y2_RK4(index_prev) 	! no vale en gral
							k2_1 = y1_RK4(index_prev) 	! no vale en gral
							
							h_improved = 0.5_dp*h
							
							y1_improved_k2 = y1_RK4(index_prev) + (k1_1*h_improved)
							y1_improved_k2 = y2_RK4(index_prev) + (k2_1*h_improved)
							k1_2 = -y2_improved_k2 	! no vale en gral
							k2_2 = y1_improved_k2 	! no vale en gral
							
							y1_improved_k3 = y1_RK4(index_prev) + (k1_2*h_improved)
							y2_improved_k3 = y2_RK4(index_prev) + (k2_2*h_improved)
							k1_3 = -y2_improved_k3 	! no vale en gral
							k2_3 = y1_improved_k3 	! no vale en gral
							
							y1_improved_k4 = y1_RK4(index_prev) + (k1_3*h)
							y2_improved_k4 = y2_RK4(index_prev) + (k2_3*h)
							k1_4 = -y2_improved_k4 	! no vale en gral
							k2_4 = y1_improved_k4 	! no vale en gral
							
							! y(i+1) = (y(i) + increment_function)
							y1_RK4(i) = y1_RK4(index_prev) + (1._dp/6._dp)*h*(k1_1 + 2._dp*(k1_2+k1_3)+k1_4)
							y2_RK4(i) = y2_RK4(index_prev) + (1._dp/6._dp)*h*(k2_1 + 2._dp*(k2_2+k2_3)+k2_4)
						end do
					! agregar casos si se quieren utilizar otras variantes de RK4
					!case(2)
					case default
						write(*,*) 'Invalid RK4 type'
				end select
			!agregar casos si se quiere usar module_functions_1D o module_functions_2D
			!case(2) 
			case default
				write(*,*) 'Invalid input type'
		end select
				
	end subroutine RK4
	
end module module_EDO_segundo_orden
