module module_EDO_segundo_orden_poincare

	use module_presition
	use module_pullen_edmonds
	
	implicit none
	
	contains
	
	!++++++++++++++++++++++++++++++++++++++++++++++++++
	! Subrutina de runge kutta para resolver cuatro
	!  ecuaciones diferenciale de primer orden.
	!++++++++++++++++++++++++++++++++++++++++++++++++++

	subroutine RK4_four_eq_poincare(n,np,a1,b1,y1_0,y2_0,y3_0,y4_0,y1_RK4_map,y3_RK4_map)
		! Data dictionary: declare calling parameter types & definitions
		integer(sp), 			intent(in)		:: n,np 				! points
		real(dp), 				intent(in)		:: y1_0 			! initial condition
		real(dp), 				intent(in)		:: y2_0 			! initial condition
		real(dp), 				intent(in)		:: y3_0 			! initial condition
		real(dp), 				intent(in)		:: y4_0 			! initial condition
		real(dp), 				intent(in)		:: a1,b1			! intervalo de validez para la solución aproximada
		real(dp), dimension(np),  intent(out)  	:: y1_RK4_map,y3_RK4_map 
		
		! Data dictionary: declare local variables types & definitions
		integer(sp) 			:: i,index_prev,counter			! index loop
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
		real(dp), dimension(n) 	:: y1_RK4 			! vector con solución aproximada
		real(dp), dimension(n)	:: y2_RK4 			! vector con solución aproximada
		real(dp), dimension(n) 	:: y3_RK4 			! vector con solución aproximada
		real(dp), dimension(n) 	:: y4_RK4 			! vector con solución aproximada
		
		! remember --> a=alpha=m2/m1; b=beta=l2/l1; c=gamma=g/l1

		h = abs(b1 - a1) * ( 1._dp / (real(n,dp)-1._dp) )
		x(1) = a1
		y1_RK4(1) = y1_0 ! initial condition of generalized coordinate 1
		y2_RK4(1) = y2_0 ! initial condition of generalized coordinate 2
		y3_RK4(1) = y3_0 ! initial condition of generalized velocity 1
		y4_RK4(1) = y4_0 ! initial condition of generalized velocity 2
	
		counter = 1_sp
		

		bucle_01: do i = 2, n, 1
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

			if (i > int(n*0.5_dp)) then
                if (y2_RK4(i)*y2_RK4(index_prev) < 0._dp) then
                    y1_RK4_map(counter) = (y1_RK4(index_prev)+y1_RK4(i))*0.5_dp
                    y3_RK4_map(counter) = (y3_RK4(index_prev)+y3_RK4(i))*0.5_dp
                    counter = counter + 1_sp
                    if (counter == np + 1_sp) then
                        exit bucle_01
                    endif
                endif 
            endif
		end do bucle_01
	end subroutine RK4_four_eq_poincare
	
end module module_EDO_segundo_orden_poincare
