module module_functions_1D
	
	use module_presition

	implicit none
	
	contains
	
	function f_1D(x,function_type)
		! Data dictionary: declare calling parameter types & definitions
		real(dp), intent(in) 	:: x 				! function variable
		integer(sp), intent(in) :: function_type 	! function_type
		real(dp) 				:: f_1D 				! function expresion
		
		select case(function_type)
			case(1)
				f_1D = exp(x) 	! increasing exponential
			case(2)
				f_1D = exp(-x) ! decreasing exponential
			case(3)
				f_1D = exp(-(x*x*0.5_dp))
			case default
				write(*,*) 'Invalid function type'
		end select
	end function f_1D
	
end module module_functions_1D
