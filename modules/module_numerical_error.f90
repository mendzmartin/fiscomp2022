module module_numerical_error
	use module_presition
	
	implicit none
	
	contains
	subroutine errors_num_integrals(I_exact, I_aprox, error, error_type)
		
		real(dp), intent(in) :: I_exact, I_aprox
		integer(sp), intent(in) :: error_type
		real(dp), intent(out) :: error
		
		select case(error_type)
		
		case (1) ! absolute error
			error = abs(I_exact - I_aprox)
		case (2) ! relative error
			error = abs((I_exact - I_aprox) * ( 1._dp / I_exact))
		case default
			write(*,*) 'Invalid error type'
		end select
		
	end subroutine errors_num_integrals

end module module_numerical_error
