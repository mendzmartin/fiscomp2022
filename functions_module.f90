module functions_module
	
	use presition_module

	implicit none

	contains
	
	function exp_func(x)
		
		! Data dictionary: declare calling parameter types & definitions
		real(dp), intent(in) 	:: x 		! function variable
		real(dp) 				:: exp_func ! function expresion
		
		exp_func = exp(x)

	end function exp_func

end module functions_module
