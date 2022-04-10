module module_double_pendulum
	use module_presition
	
	implicit none
	
	contains
	
	subroutine lagrangian_dble_pendulum(q1,q2,dq1,dq2,a,b,c,w1,w2,dw1,dw2)
	
		real(dp), intent(in) :: q1,q2,dq1,dq2,a,b,c
		real(dp), intent(out) :: w1,w2,dw1,dw2
		
		real(dp) :: dq12, dq22, q12, denom
		
		dq12=dq1*dq1
		dq22=dq2*dq2
		q12=(q1-q2)
		denom=(1._dp+a*sin(q12)*sin(q12))
		
		dw1 = (-(1._dp+a)*c*sin(q1)-a*b*dq22*sin(q12)-a*cos(q12)*(dq12*sin(q12)-c*sin(q2)))*(1._dp/denom)
		dw2 = ((1._dp+a)*(dq12*sin(q12)-c*sin(q12))+cos(q12)*((1+a)*c*sin(q1)+a*b*dq22*sin(q12)))*(1._dp/(b*denom))
		w1 = dq1
		w2 = dq2
	
	end subroutine lagrangian_dble_pendulum
	

end module module_double_pendulum
