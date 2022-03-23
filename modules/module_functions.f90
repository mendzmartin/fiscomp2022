module module_functions
	
	use module_presition

	implicit none

	contains
	function f(x,function_type)
		! Data dictionary: declare calling parameter types & definitions
		real(dp), intent(in) 	:: x 				! function variable
		integer(sp), intent(in) :: function_type 	! function_type
		real(dp) 				:: f 				! function expresion
		
		select case(function_type)
			case(1)
				f = exp(x) 	! increasing exponential
			case(2)
				f = exp(-x) ! decreasing exponential
			case default
				write(*,*) 'Invalid function type'
		end select
	end function f
	
end module module_functions
