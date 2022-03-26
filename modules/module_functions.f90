module module_functions_1var
	
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
			case(3)
				f = exp(-(x*x*0.5_dp))
			case default
				write(*,*) 'Invalid function type'
		end select
	end function f
		
		! function of 2 variables
		case(2)
			function f(x,y,function_type)
				! Data dictionary: declare calling parameter types & definitions
				real(dp), intent(in) 	:: x, y 			! function variables
				integer(sp), intent(in) :: function_type 	! function_type
				real(dp) 				:: f 				! function expresion
				
				select case(function_type)
					case(1)
						f = -x*y
					case default
						write(*,*) 'Invalid function type'
				end select
			end function f
			
		case default
			write(*,*) 'Invalid dimension type'
		end select
	end subroutine dimension_type_select
	
end module module_functions
