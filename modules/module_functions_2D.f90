module module_functions_2D
	
	use module_presition

	implicit none
	
	contains
		
	! function of 2 variables

	function f_2D(x,y,function_type)
		! Data dictionary: declare calling parameter types & definitions
		real(dp), intent(in) 	:: x, y 			! function variables
		integer(sp), intent(in) :: function_type 	! function_type
		real(dp) 				:: f_2D 			! function expresion
		
		select case(function_type)
			case(1)
				f_2D = -x*y
			case default
				write(*,*) 'Invalid function type'
		end select
	end function f_2D
	
	function f_2D_pol(x1,x2,p1,p2,a1,a2,a3,function_type)
		! Data dictionary: declare calling parameter types & definitions
		real(dp), intent(in) 	:: x1, x2 			! function variables
		real(dp), intent(in) 	:: a1,a2,a3 		! factor
		real(dp), intent(in)	:: p1,p2 			! pow
		real(dp) 				:: f_2D_pol 		! function expresion
		integer(sp), intent(in) :: function_type 	! function_type
		
		select case(function_type)
			case(1)
				f_2D_pol = a1*(x1**p1) + a2*(x2**p2) + a3
			case default
				write(*,*) 'Invalid function type'
		end select
	end function f_2D_pol
	
	subroutine lagrangian_dble_pendulum(q1,q2,dq1,dq2,a,b,c,w1,w2,dw1,dw2,function_type)
		integer(sp), intent(in) :: function_type
		real(dp), intent(in) :: q1,q2,dq1,dq2,a,b,c
		real(dp), intent(out) :: w1,w2,dw1,dw2
		real(dp) :: dq12, dq22, q12, denom
	
		select case(function_type)
				case(1)
					dq12=dq1*dq1
					dq22=dq2*dq2
					q12=(q1-q2)
					denom=(1._dp+a*sin(q12)*sin(q12))
					
					dw1 = (-(1._dp+a)*c*sin(q1)-a*b*dq22*sin(q12)-a*cos(q12)*(dq12*sin(q12)-c*sin(q2)))*(1._dp/denom)
					dw2 = (((1._dp+a)*(dq12*sin(q12)-c*sin(q2)))+(cos(q12)*((1+a)*c*sin(q1)+a*b*dq22*sin(q12))))*(1._dp/(b*denom))
					w1 = dq1
					w2 = dq2
				case default
					write(*,*) 'Invalid function type'
		end select
		
	end subroutine lagrangian_dble_pendulum
	
	
end module module_functions_2D
