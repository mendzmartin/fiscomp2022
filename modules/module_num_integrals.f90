module module_num_integrals
	use module_presition
	use module_functions_1D
	use module_gauss
	
	implicit none
	
	contains
	
	!----------------------------------------------------------
	! PURPOSE
	!----------------------------------------------------------
	! Subrutina para calcular la integral usando el método
	!   Trapezoidal
	!----------------------------------------------------------
	subroutine trapez_integ ( m, a, b, h, trapez_num_integ, function_type )
		
		! Data dictionary: declare calling parameter types & definitions
		integer(dp), intent(in) 				:: m					! cantidad puntos => m = n + 1
		integer(sp), intent(in) 				:: function_type 		! tipo de función a integrar
		real(dp), intent(in) 					:: a,b					! límites de integración
		real(dp), intent(in) 					:: h					! paso de integración
		real(dp), dimension(1,1), intent(out) 	:: trapez_num_integ		! numerical integration with trapezoidal method
		
		! Data dictionary: declare local variables types & definitions
		integer(sp) 				:: i					! index loop
		real(dp), dimension (1,m) 	:: function_vector		! column vector
		real(dp), dimension (m,1) 	:: coeff_vector			! row vector
		real(dp) 					:: x_current, x_next	! function valuation points
		
		! Execution section
		x_current = a
		
		do i = 2, m-1, 1
			coeff_vector(i,1) = 2._dp
			x_current = x_current + h
			function_vector(1,i) = f_1D(x_current,function_type)
		end do
		
		coeff_vector(1,1) = 1._dp
		coeff_vector(m,1) = 1._dp
		function_vector(1,1) = f_1D(a,function_type)
		function_vector(1,m) = f_1D(b,function_type)
		
		trapez_num_integ = h*matmul(function_vector,coeff_vector)*0.5_dp
		
	end subroutine trapez_integ
	
	!----------------------------------------------------------
	! PURPOSE
	!----------------------------------------------------------
	! Subrutina para calcular la integral usando el método de
	!   Simpson 1/3 y Simpson 3/8
	!----------------------------------------------------------
	subroutine simpson_13_integ ( m, a, b, h, simps_13_num_integ, function_type, check_value )
		
		! Data dictionary: declare calling parameter types & definitions
		integer(dp), intent(in) 				:: m					! cantidad puntos => m = n + 1
		integer(sp), intent(in) 				:: function_type 		! tipo de función a integrar
		real(dp), intent(in) 					:: a,b					! límites de integración
		real(dp), intent(in) 					:: h					! paso de integración
		real(dp), dimension(1,1), intent(out) 	:: simps_13_num_integ	! numerical integration with simpson 1/3 method
		real(dp), dimension(1,1), intent(out)	:: check_value
		
		! Data dictionary: declare local variables types & definitions
		integer(sp) 				:: i					! index loop
		real(dp), dimension (1,m) 	:: function_vector		! column vector
		real(dp), dimension (m,1) 	:: coeff_vector			! row vector
		real(dp) 					:: x_current, x_next	! function valuation points
		
		! Execution section
		x_current = a
		
		check_value = 0_dp
		
		do i = 2, m-1, 1
		
			x_current = x_current + h
			
			if ( mod(i,2) == 0 ) then
				function_vector(1,i) = f_1D(x_current,function_type)
				coeff_vector(i,1) = 4._dp
				check_value = check_value + coeff_vector(i,1)
			else
				function_vector(1,i) = f_1D(x_current,function_type)
				coeff_vector(i,1) = 2._dp
				check_value = check_value + coeff_vector(i,1)
			end if
		end do
		
		coeff_vector(1,1) = 1._dp
		check_value = check_value + coeff_vector(1,1)
		coeff_vector(m,1) = 1._dp
		check_value = check_value + coeff_vector(m,1)
		check_value = (check_value*(1._dp/3._dp) - (m-1_dp))*h 	! must be zero if the integration was ok
		
		function_vector(1,1) = f_1D(a,function_type)
		function_vector(1,m) = f_1D(b,function_type)
		
		simps_13_num_integ = h*matmul(function_vector,coeff_vector)*(1._dp/3._dp)

	end subroutine simpson_13_integ
	
	subroutine simpson_38_integ ( m, a, b, h, simps_38_num_integ, function_type )
	
	! Data dictionary: declare calling parameter types & definitions
	integer(dp), intent(in) 				:: m					! cantidad puntos => m = n + 1
	integer(sp), intent(in) 				:: function_type 		! tipo de función a integrar
	real(dp), intent(in) 					:: a,b					! límites de integración
	real(dp), intent(in) 					:: h					! paso de integración
	real(dp), dimension(1,1), intent(out) 	:: simps_38_num_integ	! numerical integration with simpson 1/3 method
	
	! Data dictionary: declare local variables types & definitions
	integer(sp) 				:: i					! index loop
	real(dp), dimension (1,m) 	:: function_vector		! column vector
	real(dp), dimension (m,1) 	:: coeff_vector			! row vector
	real(dp) 					:: x_current, x_next	! function valuation points
	
	! Execution section
	x_current = a
	
	do i = 2, m-1, 1
	
		x_current = x_current + h
	
		if ( mod(i,3) == 0 ) then
			coeff_vector(i,1) = 2._dp
			function_vector(1,i) = f_1D(x_current,function_type)
		else
			coeff_vector(i,1) = 3._dp
			function_vector(1,i) = f_1D(x_current,function_type)
		end if
	end do
	
	coeff_vector(1,1) = 1._dp
	coeff_vector(m,1) = 1._dp
	function_vector(1,1) = f_1D(a,function_type)
	function_vector(1,m) = f_1D(b,function_type)
	
	simps_38_num_integ = h*matmul(function_vector,coeff_vector)*(3._dp/8._dp)
	
	end subroutine simpson_38_integ
	
	!----------------------------------------------------------
	! PURPOSE
	!----------------------------------------------------------
	! Subrutina para calcular la integral usando el método de
	!   cuadratura de Gauss-Legendre
	!----------------------------------------------------------
	subroutine gauss_integ (m, job, a, b, gauss_num_integ, function_type)
		
		! Data dictionary: declare calling parameter types & definitions
		integer(sp), intent(in) 				:: m				! cantidad puntos => m = n + 1
		integer(sp), intent(in) 				:: function_type 	! tipo de función a integrar
		integer(sp), intent(in) 				:: job 				! type of rescale the grid points
		real(dp), intent(in) 					:: a,b				! límites de integración
		real(dp), dimension(1,1), intent(out) 	:: gauss_num_integ	! numerical integration with simpson 1/3 method
		
		! Data dictionary: declare local variables types & definitions
		integer(sp) 				:: i					! index loop
		real(dp), dimension (m,1) 	:: function_vector		! raw vector
		real(dp), dimension (m,1)	:: grid_ponit_vector	! raw vector
		real(dp), dimension (m,1) 	:: coeff_vector			! raw vector
		
		! Execution section
		call gauss(m, 0, a, b, grid_ponit_vector, coeff_vector)
		
		do i = 1, m, 1
			function_vector(i,1) = f_1D(grid_ponit_vector(i,1),function_type)
		end do
		
		gauss_num_integ = matmul(transpose(function_vector),coeff_vector)

	end subroutine gauss_integ

end module module_num_integrals
