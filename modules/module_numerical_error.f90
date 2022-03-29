!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
! Modulo para computar errores numéricos.
!
!----------------------------------------------------------

module module_numerical_error
	use module_presition
	
	implicit none
	
	contains
	
!----------------------------------------------------------
! subroutine basic_errors_num sirve para computar errores
!   básicos, por ejemplo, errores absolutos y errores relativos
!   usados en integraciones numéricas 1D y EDO de 1er orden
!----------------------------------------------------------
	subroutine basic_errors_num(val_exact, val_aprox, error, error_type)
		
		real(dp), intent(in) :: val_exact, val_aprox
		integer(sp), intent(in) :: error_type
		real(dp), intent(out) :: error
		
		select case(error_type)
		
		case (1) ! absolute error
			error = abs(val_exact - val_aprox)
		case (2) ! relative error
			error = abs((val_exact - val_aprox) * ( 1._dp / val_exact))
		case default
			write(*,*) 'Invalid error type'
		end select
		
	end subroutine basic_errors_num

end module module_numerical_error

!----------------------------------------------------------
! REFERENCES
!----------------------------------------------------------
!
!----------------------------------------------------------

!----------------------------------------------------------
! COMPILATION RUL
!----------------------------------------------------------
! gfortran -g -fcheck=all -Wall -c module_numerical_error.f90
