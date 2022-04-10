module module_double_pendulum
	use module_presition
	
	implicit none
	
	contains
	
	subroutine lagrangian_dble_pendulum(q1,q2,dq1,dq2,a,b,c,w1,w2,dw1,dw2,function_type)

		implicit none

		integer(sp), 	intent(in) 	:: function_type
		real(dp), 		intent(in) 	:: q1,q2			! generalized coordinates
		real(dp), 		intent(in) 	:: dq1,dq2			! generalized velocities
		real(dp), 		intent(in) 	:: a,b,c 			! parameters
		real(dp), 		intent(out) :: w1,w2,dw1,dw2	! first order ODEs

		real(dp) :: dq12, dq22, q12, denom, num1, num2 ! math expresions
	
		select case(function_type)
				case(1)
					dq12 = dq1*dq1
					dq22 = dq2*dq2
					q12 = (q1-q2)
					denom = (1._dp+a*sin(q12)*sin(q12))
					num1 = (-(1._dp+a)*c*sin(q1)-a*b*dq22*sin(q12)-a*cos(q12)*(dq12*sin(q12)-c*sin(q2)))
					num2 = ((1._dp+a)*(dq12*sin(q12)-c*sin(q2))+cos(q12)*((1._dp+a)*c*sin(q1)+a*b*dq22*sin(q12)))

					dw1 = num1*(1._dp/denom)
					dw2 = num2*(1._dp/(b*denom))
					w1 = dq1
					w2 = dq2
				case default
					write(*,*) 'Invalid function type'
		end select
		
	end subroutine lagrangian_dble_pendulum

	!
	! ojo! esta energía está en unidades de m1*(l1^2)
	!   donde m1 y l1 son la masa y el largo del
	!   primer pendulo más cercano al pivote
	!
	subroutine energy_dble_pendulum (theta1, theta2, w1,w2, a, b, c, T, U)

		implicit none

		real(dp), intent(in) 	:: a,b,c 			! parameters
		real(dp), intent(in) 	:: theta1, theta2 	! angles
		real(dp), intent(in) 	:: w1,w2   			! angular frecuencies
		real(dp), intent(out) 	:: T, U 			! kinetic and potential energies

		real(dp) :: w12,w22,b2,theta12 ! math expresions	

		w12 = w1*w1
		w22	= w2*w2
		theta12 = (theta1 - theta2)
		b2	= b*b

		T = 0.5_dp*((1._dp+a)*w12+a*b2*w22+2._dp*a*b*w1*w2*cos(theta12))
		U = -c*((1._dp+a)*cos(theta1)+a*b*cos(theta2))

	end subroutine energy_dble_pendulum

end module module_double_pendulum
