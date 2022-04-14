module module_EDO_segundo_orden

	use module_presition
	use module_functions_2D
	use module_double_pendulum
	use module_pullen_edmonds
	
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
		case(1) ! using f_2D_pol(x1,x2,p1,p2,a1,a2,a3,function_type)
			do i = 2, n, 1
				index_prev = i-1
				x(i) = x(index_prev) + h
				y1_euler(i) = y1_euler(index_prev) + h*f_2D_pol(0._dp,y2_euler(index_prev),0._dp,1._dp,0._dp,-1._dp,0._dp,1)
				y2_euler(i) = y2_euler(index_prev) + h*f_2D_pol(0._dp,y1_euler(index_prev),0._dp,1._dp,0._dp,1._dp,0._dp,1)
			end do
		!case(2) ! add to use f_2D(x,y,function_type)
		!case(3) ! add to use f_1D(x,function_type)
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
			case(1) ! using f_2D_pol(x,y,p1,p2,a1,a2,a3,function_type)
				select case (RK2_type)
					case(1) ! Método de Heun con un solo corrector (a_2 = 1/2)
						h = abs(b - a) * ( 1._dp / (n-1_sp) )
						x(1) = a
						y1_RK2(1) = y1_0
						y2_RK2(1) = y2_0
					
						do i = 2, n, 1
							index_prev = i-1
							x(i) = x(index_prev) + h
							
							! ki_1 = fi(y)
							k1_1 = f_2D_pol(0._dp,y2_RK2(index_prev),0._dp,1._dp,0._dp,-1._dp,0._dp,1)
							k2_1 = f_2D_pol(0._dp,y1_RK2(index_prev),0._dp,1._dp,0._dp,1._dp,0._dp,1)
							
							! yi_improved = yi_RK2(index_prev) + ki_1*h
							y1_improved = y1_RK2(index_prev) + (k1_1*h)
							y2_improved = y2_RK2(index_prev) + (k2_1*h)
							
							! ki_2 = fi(y_improve)
							k1_2 = f_2D_pol(0._dp,y2_improved,0._dp,1._dp,0._dp,-1._dp,0._dp,1)
							k2_2 = f_2D_pol(0._dp,y1_improved,0._dp,1._dp,0._dp,1._dp,0._dp,1)
							
							! y(i+1) = (y(i) + increment_function)
							y1_RK2(i) = y1_RK2(index_prev) + (0.5_dp)*h*(k1_1 + k1_2)
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
							
							k1_1 = f_2D_pol(0._dp,y2_RK2(index_prev),0._dp,1._dp,0._dp,-1._dp,0._dp,1)
							k2_1 = f_2D_pol(0._dp,y1_RK2(index_prev),0._dp,1._dp,0._dp,1._dp,0._dp,1)
							h_improved = 0.5_dp*h
							y1_improved = y1_RK2(index_prev) + (k1_1*h_improved)
							y2_improved = y2_RK2(index_prev) + (k2_1*h_improved)
							k1_2 = f_2D_pol(0._dp,y2_improved,0._dp,1._dp,0._dp,-1._dp,0._dp,1)
							k2_2 = f_2D_pol(0._dp,y1_improved,0._dp,1._dp,0._dp,1._dp,0._dp,1)
							
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
							
							k1_1 = f_2D_pol(0._dp,y2_RK2(index_prev),0._dp,1._dp,0._dp,-1._dp,0._dp,1)
							k2_1 = f_2D_pol(0._dp,y1_RK2(index_prev),0._dp,1._dp,0._dp,1._dp,0._dp,1)
							h_improved = 0.25_dp*h
							y1_improved = y1_RK2(index_prev) + (3._dp*k1_1*h_improved)
							y2_improved = y2_RK2(index_prev) + (3._dp*k2_1*h_improved)
							k1_2 = f_2D_pol(0._dp,y2_improved,0._dp,1._dp,0._dp,-1._dp,0._dp,1)
							k2_2 = f_2D_pol(0._dp,y1_improved,0._dp,1._dp,0._dp,1._dp,0._dp,1)
							
							! y(i+1) = (y(i) + increment_function)
							y1_RK2(i) = y1_RK2(index_prev) + (1._dp/3._dp)*h*(k1_1+2*k1_2)
							y2_RK2(i) = y2_RK2(index_prev) + (1._dp/3._dp)*h*(k2_1+2*k2_2)
						end do
					case default
						write(*,*) 'Invalid RK2 type'
				end select
			!case(2) ! add to use f_2D(x,y,function_type)
			!case(3) ! add to use f_1D(x,function_type)
			case default
				write(*,*) 'Invalid input type'
		end select
	end subroutine RK2
	
	!++++++++++++++++++++++++++++++++++++++++++++++++++
	! Subrutina de runge kutta para resolver dos
	!  ecuaciones diferenciale de primer orden.
	!++++++++++++++++++++++++++++++++++++++++++++++++++

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
			case(1) ! using f_2D_pol(x1,x2,p1,p2,a1,a2,a3,function_type)
				select case (RK4_type)
					case(1) ! Método clásico
						h = abs(b - a) * ( 1._dp / (n-1_sp) )
						x(1) = a
						y1_RK4(1) = y1_0
						y2_RK4(1) = y2_0
					
						do i = 2, n, 1
							index_prev = i-1
							x(i) = x(index_prev) + h
							
							k1_1 = f_2D_pol(0._dp,y2_RK4(index_prev),0._dp,1._dp,0._dp,-1._dp,0._dp,1)
							k2_1 = f_2D_pol(0._dp,y1_RK4(index_prev),0._dp,1._dp,0._dp,1._dp,0._dp,1)
							
							h_improved = 0.5_dp*h
							
							y1_improved_k2 = y1_RK4(index_prev) + (k1_1*h_improved)
							y2_improved_k2 = y2_RK4(index_prev) + (k2_1*h_improved)
							k1_2 = f_2D_pol(0._dp,y2_improved_k2,0._dp,1._dp,0._dp,-1._dp,0._dp,1)
							k2_2 = f_2D_pol(0._dp,y1_improved_k2,0._dp,1._dp,0._dp,1._dp,0._dp,1)
							
							y1_improved_k3 = y1_RK4(index_prev) + (k1_2*h_improved)
							y2_improved_k3 = y2_RK4(index_prev) + (k2_2*h_improved)
							k1_3 = f_2D_pol(0._dp,y2_improved_k3,0._dp,1._dp,0._dp,-1._dp,0._dp,1)
							k2_3 = f_2D_pol(0._dp,y1_improved_k3,0._dp,1._dp,0._dp,1._dp,0._dp,1)
							
							y1_improved_k4 = y1_RK4(index_prev) + (k1_3*h)
							y2_improved_k4 = y2_RK4(index_prev) + (k2_3*h)
							k1_4 = f_2D_pol(0._dp,y2_improved_k4,0._dp,1._dp,0._dp,-1._dp,0._dp,1)
							k2_4 = f_2D_pol(0._dp,y1_improved_k4,0._dp,1._dp,0._dp,1._dp,0._dp,1)
							
							! y(i+1) = (y(i) + increment_function)
							y1_RK4(i) = y1_RK4(index_prev) + (1._dp/6._dp)*h*(k1_1 + 2._dp*(k1_2+k1_3)+k1_4)
							y2_RK4(i) = y2_RK4(index_prev) + (1._dp/6._dp)*h*(k2_1 + 2._dp*(k2_2+k2_3)+k2_4)
						end do
					case default
						write(*,*) 'Invalid RK4 type'
				end select
			!case(2) ! add to use f_2D(x,y,function_type)
			!case(3) ! add to use f_1D(x,function_type)
			case default
				write(*,*) 'Invalid input type'
		end select
				
	end subroutine RK4
	
	!++++++++++++++++++++++++++++++++++++++++++++++++++
	! Subrutina de runge kutta para resolver cuatro
	!  ecuaciones diferenciale de primer orden.
	!++++++++++++++++++++++++++++++++++++++++++++++++++

	subroutine RK4_four_eq(RK4_type,n,a1,b1,y1_0,y2_0,y3_0,y4_0,y1_RK4,y2_RK4,y3_RK4,y4_RK4,function_type,input_type)
		! Data dictionary: declare calling parameter types & definitions
		integer(sp), 			intent(in)		:: RK4_type
		integer(sp), 			intent(in)		:: n 				! points
		integer(sp), 			intent(in) 		:: function_type 	! tipo de función a integrar
		integer(sp), 			intent(in)		:: input_type 		! type of input to create function
		real(dp), 				intent(in)		:: y1_0 			! initial condition
		real(dp), 				intent(in)		:: y2_0 			! initial condition
		real(dp), 				intent(in)		:: y3_0 			! initial condition
		real(dp), 				intent(in)		:: y4_0 			! initial condition
		real(dp), 				intent(in)		:: a1,b1			! intervalo de validez para la solución aproximada
		real(dp), dimension(n), intent(out) 	:: y1_RK4 			! vector con solución aproximada
		real(dp), dimension(n), intent(out) 	:: y2_RK4 			! vector con solución aproximada
		real(dp), dimension(n), intent(out) 	:: y3_RK4 			! vector con solución aproximada
		real(dp), dimension(n), intent(out) 	:: y4_RK4 			! vector con solución aproximada
		
		! Data dictionary: declare local variables types & definitions
		integer(sp) 			:: i,index_prev 			! index loop
		real(dp) 				:: k1_1, k1_2, k1_3, k1_4	! recurrence relations
		real(dp) 				:: k2_1, k2_2, k2_3, k2_4	! recurrence relations
		real(dp) 				:: k3_1, k3_2, k3_3, k3_4	! recurrence relations
		real(dp) 				:: k4_1, k4_2, k4_3, k4_4	! recurrence relations
		real(dp) 				:: y1_improved_k2			! improved function evaluation
		real(dp) 				:: y1_improved_k3			! improved function evaluation
		real(dp) 				:: y1_improved_k4			! improved function evaluation
		real(dp) 				:: y2_improved_k2			! improved function evaluation
		real(dp) 				:: y2_improved_k3			! improved function evaluation
		real(dp) 				:: y2_improved_k4			! improved function evaluation
		real(dp) 				:: y3_improved_k2			! improved function evaluation
		real(dp) 				:: y3_improved_k3			! improved function evaluation
		real(dp) 				:: y3_improved_k4			! improved function evaluation
		real(dp) 				:: y4_improved_k2			! improved function evaluation
		real(dp) 				:: y4_improved_k3			! improved function evaluation
		real(dp) 				:: y4_improved_k4			! improved function evaluation
		real(dp) 				:: h_improved 				! improved step evaluation
		real(dp), dimension(n) 	:: x 						! grid points vector
		real(dp) 				:: h 						! step
		real(dp) 				:: f1,f2,f3,f4 				! custom functions
		
		select case(input_type)
			case(1) ! using lagrangian_dble_pendulum(q1,q2,dq1,dq2,a,b,c,w1,w2,dw1,dw2,function_type)
				select case (RK4_type)
					case(1) ! Método clásico

						! remember --> a=alpha=m2/m1; b=beta=l2/l1; c=gamma=g/l1

						h = abs(b1 - a1) * ( 1._dp / (real(n,dp)-1._dp) )
						x(1) = a1
						y1_RK4(1) = y1_0 ! initial condition of generalized coordinate 1
						y2_RK4(1) = y2_0 ! initial condition of generalized coordinate 2
						y3_RK4(1) = y3_0 ! initial condition of generalized velocity 1
						y4_RK4(1) = y4_0 ! initial condition of generalized velocity 2
					
						do i = 2, n, 1
							index_prev = i-1
							x(i) = x(index_prev) + h
					
							call lagrangian_dble_pendulum(y1_RK4(index_prev),y2_RK4(index_prev),y3_RK4(index_prev),&
							y4_RK4(index_prev),(1._dp/3._dp),0.5_dp,0.5_dp,f1,f2,f3,f4,1_sp)
					
							k1_1 = f1
							k2_1 = f2
							k3_1 = f3
							k4_1 = f4
							
							h_improved = 0.5_dp*h
							
							y1_improved_k2 = y1_RK4(index_prev) + (k1_1*h_improved)
							y2_improved_k2 = y2_RK4(index_prev) + (k2_1*h_improved)
							y3_improved_k2 = y3_RK4(index_prev) + (k3_1*h_improved)
							y4_improved_k2 = y4_RK4(index_prev) + (k4_1*h_improved)
							
							call lagrangian_dble_pendulum(y1_improved_k2,y2_improved_k2,y3_improved_k2,&
							y4_improved_k2,(1._dp/3._dp),0.5_dp,0.5_dp,f1,f2,f3,f4,1_sp)
							
							k1_2 = f1
							k2_2 = f2
							k3_2 = f3
							k4_2 = f4
							
							y1_improved_k3 = y1_RK4(index_prev) + (k1_2*h_improved)
							y2_improved_k3 = y2_RK4(index_prev) + (k2_2*h_improved)
							y3_improved_k3 = y3_RK4(index_prev) + (k3_2*h_improved)
							y4_improved_k3 = y4_RK4(index_prev) + (k4_2*h_improved)
							
							call lagrangian_dble_pendulum(y1_improved_k3,y2_improved_k3,y3_improved_k3,&
							y4_improved_k3,(1._dp/3._dp),0.5_dp,0.5_dp,f1,f2,f3,f4,1_sp)
							k1_3 = f1
							k2_3 = f2
							k3_3 = f3
							k4_3 = f4
							
							y1_improved_k4 = y1_RK4(index_prev) + (k1_3*h)
							y2_improved_k4 = y2_RK4(index_prev) + (k2_3*h)
							y3_improved_k4 = y3_RK4(index_prev) + (k3_3*h)
							y4_improved_k4 = y4_RK4(index_prev) + (k4_3*h)
							
							call lagrangian_dble_pendulum(y1_improved_k4,y2_improved_k4,y3_improved_k4,&
							y4_improved_k4,(1._dp/3._dp),0.5_dp,0.5_dp,f1,f2,f3,f4,1_sp)
							k1_4 = f1
							k2_4 = f2
							k3_4 = f3
							k4_4 = f4
							
							! y(i+1) = (y(i) + increment_function)
							y1_RK4(i) = y1_RK4(index_prev) + (1._dp/6._dp)*h*(k1_1 + 2._dp*(k1_2+k1_3)+k1_4)
							y2_RK4(i) = y2_RK4(index_prev) + (1._dp/6._dp)*h*(k2_1 + 2._dp*(k2_2+k2_3)+k2_4)
							y3_RK4(i) = y3_RK4(index_prev) + (1._dp/6._dp)*h*(k3_1 + 2._dp*(k3_2+k3_3)+k3_4)
							y4_RK4(i) = y4_RK4(index_prev) + (1._dp/6._dp)*h*(k4_1 + 2._dp*(k4_2+k4_3)+k4_4)
						end do
					case default
						write(*,*) 'Invalid RK4 type'
				end select
			case(2) ! using hamiltonian_pullen_edmonds(q1,q2,p1,p2,alpha,m,omega,dq1,dq2,dp1,dp2,function_type)
				select case (RK4_type)
				case(1) ! Método clásico
					! remember --> a=alpha=m2/m1; b=beta=l2/l1; c=gamma=g/l1

					h = abs(b1 - a1) * ( 1._dp / (real(n,dp)-1._dp) )
					x(1) = a1
					y1_RK4(1) = y1_0 ! initial condition of generalized coordinate 1
					y2_RK4(1) = y2_0 ! initial condition of generalized coordinate 2
					y3_RK4(1) = y3_0 ! initial condition of generalized velocity 1
					y4_RK4(1) = y4_0 ! initial condition of generalized velocity 2
				
					do i = 2, n, 1
						index_prev = i-1
						x(i) = x(index_prev) + h
				
						call hamiltonian_pullen_edmonds(y1_RK4(index_prev),y2_RK4(index_prev),y3_RK4(index_prev),&
						y4_RK4(index_prev),0.05_dp,1._dp,1._dp,f1,f2,f3,f4,1_sp)
				
						k1_1 = f1
						k2_1 = f2
						k3_1 = f3
						k4_1 = f4
						
						h_improved = 0.5_dp*h
						
						y1_improved_k2 = y1_RK4(index_prev) + (k1_1*h_improved)
						y2_improved_k2 = y2_RK4(index_prev) + (k2_1*h_improved)
						y3_improved_k2 = y3_RK4(index_prev) + (k3_1*h_improved)
						y4_improved_k2 = y4_RK4(index_prev) + (k4_1*h_improved)
						
						call hamiltonian_pullen_edmonds(y1_improved_k2,y2_improved_k2,y3_improved_k2,&
						y4_improved_k2,0.05_dp,1._dp,1._dp,f1,f2,f3,f4,1_sp)
						
						k1_2 = f1
						k2_2 = f2
						k3_2 = f3
						k4_2 = f4
						
						y1_improved_k3 = y1_RK4(index_prev) + (k1_2*h_improved)
						y2_improved_k3 = y2_RK4(index_prev) + (k2_2*h_improved)
						y3_improved_k3 = y3_RK4(index_prev) + (k3_2*h_improved)
						y4_improved_k3 = y4_RK4(index_prev) + (k4_2*h_improved)
						
						call hamiltonian_pullen_edmonds(y1_improved_k3,y2_improved_k3,y3_improved_k3,&
						y4_improved_k3,0.05_dp,1._dp,1._dp,f1,f2,f3,f4,1_sp)

						k1_3 = f1
						k2_3 = f2
						k3_3 = f3
						k4_3 = f4
						
						y1_improved_k4 = y1_RK4(index_prev) + (k1_3*h)
						y2_improved_k4 = y2_RK4(index_prev) + (k2_3*h)
						y3_improved_k4 = y3_RK4(index_prev) + (k3_3*h)
						y4_improved_k4 = y4_RK4(index_prev) + (k4_3*h)
						
						call hamiltonian_pullen_edmonds(y1_improved_k4,y2_improved_k4,y3_improved_k4,&
						y4_improved_k4,0.05_dp,1._dp,1._dp,f1,f2,f3,f4,1_sp)

						k1_4 = f1
						k2_4 = f2
						k3_4 = f3
						k4_4 = f4
						
						! y(i+1) = (y(i) + increment_function)
						y1_RK4(i) = y1_RK4(index_prev) + (1._dp/6._dp)*h*(k1_1 + 2._dp*(k1_2+k1_3)+k1_4)
						y2_RK4(i) = y2_RK4(index_prev) + (1._dp/6._dp)*h*(k2_1 + 2._dp*(k2_2+k2_3)+k2_4)
						y3_RK4(i) = y3_RK4(index_prev) + (1._dp/6._dp)*h*(k3_1 + 2._dp*(k3_2+k3_3)+k3_4)
						y4_RK4(i) = y4_RK4(index_prev) + (1._dp/6._dp)*h*(k4_1 + 2._dp*(k4_2+k4_3)+k4_4)
					end do
				case default
					write(*,*) 'Invalid RK4 type'
				end select
			case default
				write(*,*) 'Invalid input type'
		end select	
	end subroutine RK4_four_eq
	
end module module_EDO_segundo_orden
