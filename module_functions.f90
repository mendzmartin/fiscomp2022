module module_functions
	
	use module_presition

	implicit none

	contains
	
	function exp_func(x)
		
		! Data dictionary: declare calling parameter types & definitions
		real(dp), intent(in) 	:: x 		! function variable
		real(dp) 				:: exp_func ! function expresion
		
		exp_func = exp(x)

	end function exp_func

end module module_functions
