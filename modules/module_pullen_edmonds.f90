module module_pullen_edmonds
	use module_presition
	
	implicit none
	
	contains
	
	!subroutine hamiltonian_pullen_edmonds(q1,q2,p1,p2,alpha,dq1,dq2,dp1,dp2,function_type)

	subroutine hamiltonian_pullen_edmonds(q1,q2,p1,p2,alpha,m,omega,dq1,dq2,dp1,dp2,function_type)
		integer(sp), 	intent(in) 	:: function_type
		real(dp), 		intent(in) 	:: q1,q2			! generalized coordinates
		real(dp), 		intent(in) 	:: p1,p2			! generalized velocities
		real(dp), 		intent(in) 	:: alpha 			! parameter
		real(dp), 		intent(out) :: dq1,dq2,dp1,dp2	! first order ODEs
		real(dp), 		intent(in) 	:: m,omega ! mass and angular frecuency

		real(dp) :: q12, q22, omega2
		q12=q1*q1
		q22=q2*q2
		omega2=omega*omega

		select case(function_type)
			case(1) ! funcion custom
				dq1=p1*(1._dp/m)
				dq2=p2*(1._dp/m)
				dp1=-q1*(m*omega2+2._dp*alpha*q22)
				dp2=-q2*(m*omega2+2._dp*alpha*q12)
			case default
				write(*,*) 'Invalid function type'
		end select
	end subroutine hamiltonian_pullen_edmonds

    subroutine energy_pullen_edmons(q1,q2,p1,p2,alpha,m,omega,E,function_type)
        integer(sp), 	intent(in) 	:: function_type
		real(dp), 		intent(in) 	:: q1,q2			! generalized coordinates
		real(dp), 		intent(in) 	:: p1,p2			! generalized velocities
		real(dp), 		intent(in) 	:: alpha 			! parameter
		real(dp), 		intent(out) :: E	            ! first order ODEs
		real(dp), 		intent(in) 	:: m,omega          ! mass and angular frecuency
    
		real(dp) :: p12, p22, q12, q22, omega2
        p12=p1*p1
        p22=p2*p2
		q12=q1*q1
		q22=q2*q2
		omega2=omega*omega

        E = 0.5_dp*(1._dp/m)*(p12+p22)+0.5_dp*m*omega2*(q12+q22)+alpha*q12*q22
        
    end subroutine energy_pullen_edmons

	subroutine q2_pullen_edmons(q1,q2,p1,p2,alpha,m,omega,E,control)
		real(dp), 		intent(in) 	:: E,q1			! generalized coordinates
		real(dp), 		intent(in) 	:: p1,p2			! generalized velocities
		real(dp), 		intent(in) 	:: alpha 			! parameter
		real(dp), 		intent(out) :: q2	            ! first order ODEs
		real(dp), 		intent(in) 	:: m,omega          ! mass and angular frecuency
		logical, intent(out) :: control
    
		real(dp) :: p12, p22, q12, omega2, factor
        p12=p1*p1
        p22=p2*p2
		q12=q1*q1
		omega2=omega*omega

		factor = 0.5_dp*(p12+p22)*(1._dp/m)-0.5_dp*m*omega2*q12
		if (E >= factor) then
			control = .true.
        	q2=sqrt(E-factor)*(1._dp/(m*omega2*0.5_dp+alpha*q12))
		else
			control = .false.
			q2=0
		end if
        
    end subroutine q2_pullen_edmons

end module module_pullen_edmonds
