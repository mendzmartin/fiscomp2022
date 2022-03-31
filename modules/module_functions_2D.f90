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
	
end module module_functions_2D
