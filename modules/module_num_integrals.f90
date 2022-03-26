module module_num_integrals
	use module_presition
	use module_functions
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
		
		do i = 2, m-2, 1
			coeff_vector(i,1) = 2._dp
			function_vector(1,i) = f(x_current,function_type)

			x_current = x_current + h
		end do
		
		coeff_vector(1,1) = 1._dp
		coeff_vector(m-1,1) = 1._dp
		function_vector(1,1) = f(a,function_type)
		function_vector(1,m-1) = f(b,function_type)
		
		trapez_num_integ = h*matmul(function_vector,coeff_vector)*0.5_dp
		
	end subroutine trapez_integ
	
	!----------------------------------------------------------
	! PURPOSE
	!----------------------------------------------------------
	! Subrutina para calcular la integral usando el método de
	!   Simpson 1/3 y Simpson 3/8
	!----------------------------------------------------------
	subroutine simpson_13_integ ( m, a, b, h, simps_13_num_integ, function_type )
		
		! Data dictionary: declare calling parameter types & definitions
		integer(dp), intent(in) 				:: m					! cantidad puntos => m = n + 1
		integer(sp), intent(in) 				:: function_type 		! tipo de función a integrar
		real(dp), intent(in) 					:: a,b					! límites de integración
		real(dp), intent(in) 					:: h					! paso de integración
		real(dp), dimension(1,1), intent(out) 	:: simps_13_num_integ	! numerical integration with simpson 1/3 method
		
		! Data dictionary: declare local variables types & definitions
		integer(sp) 				:: i					! index loop
		real(dp), dimension (1,m) 	:: function_vector		! column vector
		real(dp), dimension (m,1) 	:: coeff_vector			! row vector
		real(dp) 					:: x_current, x_next	! function valuation points
		
		! Execution section
		x_current = a
		
		do i = 2, m-2, 1
			if ( mod(i,2) == 0 ) then
				coeff_vector(i,1) = 4._dp
				function_vector(1,i) = f(x_current,function_type)
			else
				coeff_vector(i,1) = 2._dp
				function_vector(1,i) = f(x_current,function_type)
			end if
			
			x_current = x_current + h
		end do
		
		coeff_vector(1,1) = 1._dp
		coeff_vector(m-1,1) = 1._dp
		function_vector(1,1) = f(a,function_type)
		function_vector(1,m-1) = f(b,function_type)
		
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
	
	do i = 2, m-2, 1
		if ( mod(i,3) == 0 ) then
			coeff_vector(i,1) = 2._dp
			function_vector(1,i) = f(x_current,function_type)
		else
			coeff_vector(i,1) = 3._dp
			function_vector(1,i) = f(x_current,function_type)
		end if
		
		x_current = x_current + h
	end do
	
	coeff_vector(1,1) = 1._dp
	coeff_vector(m-1,1) = 1._dp
	function_vector(1,1) = f(a,function_type)
	function_vector(1,m-1) = f(b,function_type)
	
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
			function_vector(i,1) = f(grid_ponit_vector(i,1),function_type)
		end do
		
		gauss_num_integ = matmul(transpose(function_vector),coeff_vector)

	end subroutine gauss_integ

end module module_num_integrals
