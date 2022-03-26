module module_functions_2D
	
	use module_presition

	implicit none
	
	contains
		
	! function of 2 variables

	function f_2D(x,y,function_type)
		! Data dictionary: declare calling parameter types & definitions
		real(dp), intent(in) 	:: x, y 			! function variables
		integer(sp), intent(in) :: function_type 	! function_type
		real(dp) 				:: f_2D 				! function expresion
		
		select case(function_type)
			case(1)
				f_2D = -x*y
			case default
				write(*,*) 'Invalid function type'
		end select
	end function f_2D
	
end module module_functions_2D
