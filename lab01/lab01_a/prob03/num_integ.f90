!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
! Integración numérica.
! Puntos equiespaciados: Comparación Trapezoidal vs. Simpson
! 	programa sólo válido para la integral definida de la
! 	función exponencial f(x) = exp(x)
!----------------------------------------------------------

program num_integ

	implicit none

	! Explicit variables declaration
	integer, parameter 			:: sp = selected_real_kind( p=6, r=37 )		! simple presicion (sp) class
	integer, parameter 			:: dp = selected_real_kind( p=15, r=307 )	! double presicion (dp) class
	real(dp)					:: exact_integ								! valor de la integral exacta
	real(dp), dimension(1,1)	:: trapez_num_integ							! numerical integration with trapezoidal method
	real(dp), dimension(1,1) 	:: simps_13_num_integ						! numerical integration with simpson 1/3 method
	real(dp), dimension(1,1) 	:: simps_38_num_integ						! numerical integration with simpson 3/8 method
	integer(dp)					:: k										! potencia
	integer(dp) 				:: n, m 									! número de intervalos y número de puntos
	real(dp)					:: a,b										! límites de integración
	real(dp)					:: h										! paso de integración
	
	! Datos pedidos al usuario
	write( *, * ) 'Ingrese el limite inferior de integración (a) y presione Enter.'
	read( *, * ) a
	write( *, * ) 'Ingrese el limite superior de integración (a) y presione Enter.'
	read( *, * ) b
	write( *, * ) '---------------------------------------'
			
	20 format (A15, F15.6)
	21 format (A5, I8, A10, F8.6)
	
	do k = 2, 8, 2
	
		n = 2d0**k
		m = n + 1d0 
		h = ( b - a ) * ( 1d0 / n )
		exact_integ = ( exp( b ) - exp(a) )
		
		call trapez_integ ( m, a, b, h, trapez_num_integ )
		call simpson_13_integ ( m, a, b, h, simps_13_num_integ )
		call simpson_38_integ ( m, a, b, h, simps_38_num_integ )
		
		write( *, 21 ) 'n = ', n, 'h = ', h
		write( *, 20 ) 'I_exact       = ', exact_integ
		write( *, 20 ) 'I_trapezoidal = ', trapez_num_integ
		write( *, 20 ) 'I_simpson13   = ', simps_13_num_integ
		write( *, 20 ) 'I_simpson38   = ', simps_38_num_integ
		write( *, * ) '---------------------------------------'
	end do
	
	return
end program num_integ



!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
! Subrutina para calcular la integral usando el método
!   Trapezoidal
!----------------------------------------------------------
subroutine trapez_integ ( m, a, b, h, trapez_num_integ )

	implicit none
	
	integer, parameter 	:: sp = selected_real_kind( p=6, r=37 )		! simple presicion (sp) class
	integer, parameter	:: dp = selected_real_kind( p=15, r=307 )	! double presicion (dp) class
	
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
		function_vector(1,i) = exp(x_current)
		x_current = x_current + h
	end do
	
	coeff_vector(1,1) = 1.d0
	coeff_vector(m-1,1) = 1.d0
	function_vector(1,1) = exp(a)
	function_vector(1,m-1) = exp(b)
	
	trapez_num_integ = h*matmul(function_vector,coeff_vector)*0.5d0
	
	return
end subroutine trapez_integ

!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
! Subrutina para calcular la integral usando el método de
!   Simpson 1/3 y Simpson 3/8
!----------------------------------------------------------
subroutine simpson_13_integ ( m, a, b, h, simps_13_num_integ )
	implicit none
	
	integer, parameter 	:: sp = selected_real_kind( p=6, r=37 )		! simple presicion (sp) class
	integer, parameter	:: dp = selected_real_kind( p=15, r=307 )	! double presicion (dp) class
	
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
			function_vector(1,i) = exp(x_current)
		else
			coeff_vector(i,1) = 2d0
			function_vector(1,i) = exp(x_current)
		end if
		
		x_current = x_current + h
	end do
	
	coeff_vector(1,1) = 1.d0
	coeff_vector(m-1,1) = 1.d0
	function_vector(1,1) = exp(a)
	function_vector(1,m-1) = exp(b)
	
	simps_13_num_integ = h*matmul(function_vector,coeff_vector)*(1d0/3d0)
	
	return
end subroutine simpson_13_integ

subroutine simpson_38_integ ( m, a, b, h, simps_38_num_integ )
	implicit none
	
	integer, parameter 	:: sp = selected_real_kind( p=6, r=37 )		! simple presicion (sp) class
	integer, parameter	:: dp = selected_real_kind( p=15, r=307 )	! double presicion (dp) class
	
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
			function_vector(1,i) = exp(x_current)
		else
			coeff_vector(i,1) = 3d0
			function_vector(1,i) = exp(x_current)
		end if
		
		x_current = x_current + h
	end do
	
	coeff_vector(1,1) = 1.d0
	coeff_vector(m-1,1) = 1.d0
	function_vector(1,1) = exp(a)
	function_vector(1,m-1) = exp(b)
	
	simps_38_num_integ = h*matmul(function_vector,coeff_vector)*(3d0/8d0)
	
	return
end subroutine simpson_38_integ


!----------------------------------------------------------
! REFERENCES
!----------------------------------------------------------
!
!
!
!----------------------------------------------------------

!----------------------------------------------------------
! COMPILATION RUL
!----------------------------------------------------------
! gfortran -g3 -Wall -Wextra -Wconversion num_integ.f90 -o num_integ && ./num_integ
!----------------------------------------------------------
