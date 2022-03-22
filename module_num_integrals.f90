module module_num_integrals
	use module_presition
	use module_functions
	
	implicit none
	
	contains
	
	!----------------------------------------------------------
	! PURPOSE
	!----------------------------------------------------------
	! Subrutina para calcular la integral usando el método
	!   Trapezoidal
	!----------------------------------------------------------
	subroutine trapez_integ ( m, a, b, h, trapez_num_integ )
		
		! Data dictionary: declare calling parameter types & definitions
		integer(dp), intent(in) 				:: m					! cantidad puntos => m = n + 1
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
			coeff_vector(i,1) = 2d0
			function_vector(1,i) = exp_func(x_current)
			x_current = x_current + h
		end do
		
		coeff_vector(1,1) = 1.d0
		coeff_vector(m-1,1) = 1.d0
		function_vector(1,1) = exp_func(a)
		function_vector(1,m-1) = exp_func(b)
		
		trapez_num_integ = h*matmul(function_vector,coeff_vector)*0.5d0
		
	end subroutine trapez_integ
	
	!----------------------------------------------------------
	! PURPOSE
	!----------------------------------------------------------
	! Subrutina para calcular la integral usando el método de
	!   Simpson 1/3 y Simpson 3/8
	!----------------------------------------------------------
	subroutine simpson_13_integ ( m, a, b, h, simps_13_num_integ )
		
		! Data dictionary: declare calling parameter types & definitions
		integer(dp), intent(in) 				:: m					! cantidad puntos => m = n + 1
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
				coeff_vector(i,1) = 4d0
				function_vector(1,i) = exp_func(x_current)
			else
				coeff_vector(i,1) = 2d0
				function_vector(1,i) = exp_func(x_current)
			end if
			
			x_current = x_current + h
		end do
		
		coeff_vector(1,1) = 1.d0
		coeff_vector(m-1,1) = 1.d0
		function_vector(1,1) = exp_func(a)
		function_vector(1,m-1) = exp_func(b)
		
		simps_13_num_integ = h*matmul(function_vector,coeff_vector)*(1d0/3d0)

	end subroutine simpson_13_integ
	
	subroutine simpson_38_integ ( m, a, b, h, simps_38_num_integ )
	
	! Data dictionary: declare calling parameter types & definitions
	integer(dp), intent(in) 				:: m					! cantidad puntos => m = n + 1
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
			coeff_vector(i,1) = 2d0
			function_vector(1,i) = exp_func(x_current)
		else
			coeff_vector(i,1) = 3d0
			function_vector(1,i) = exp_func(x_current)
		end if
		
		x_current = x_current + h
	end do
	
	coeff_vector(1,1) = 1.d0
	coeff_vector(m-1,1) = 1.d0
	function_vector(1,1) = exp_func(a)
	function_vector(1,m-1) = exp_func(b)
	
	simps_38_num_integ = h*matmul(function_vector,coeff_vector)*(3d0/8d0)
	
	end subroutine simpson_38_integ

end module module_num_integrals
