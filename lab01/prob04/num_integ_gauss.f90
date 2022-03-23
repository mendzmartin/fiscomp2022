!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
! Integración numérica.
! Puntos equiespaciados: Comparación Trapezoidal vs. Simpson
! 	programa sólo válido para la integral definida de la
!----------------------------------------------------------

program num_integ_gauss
	
	use module_presition 							! use module of presition
	use module_functions 							! use module of functions
	use module_num_integrals 						! use module of numerical integrals
	use module_gauss 								! use module of gauss-legendre cuadrature
	use module_numerical_error 	! use module of numerical errors
	
	implicit none

	! Explicit variables declaration
	real(dp)					:: exact_integ			! valor de la integral exacta
	real(dp), dimension(1,1)	:: trapez_num_integ		! numerical integration with trapezoidal method
	real(dp), dimension(1,1) 	:: simps_13_num_integ	! numerical integration with simpson 1/3 method
	real(dp), dimension(1,1) 	:: simps_38_num_integ	! numerical integration with simpson 3/8 method
	real(dp), dimension(1,1) 	:: gauss_num_integ		! numerical integration with gauss quadrature method
	integer(dp)					:: k					! potencia
	integer(dp) 				:: n, m 				! número de intervalos y número de puntos
	real(dp)					:: a,b					! límites de integración
	real(dp)					:: h					! paso de integración
	integer(sp) 				:: m_gauss 				! número de intervalos para usar en modulo de gauss
	integer(sp)					:: function_type 		! tipo de función a integrar
	real(dp) 					:: rel_error 			! error relativo
	
	! Datos pedidos al usuario
	write( *, * ) 'Ingrese el limite inferior de integración (a) y presione Enter.'
	read( *, * ) a
	write( *, * ) 'Ingrese el limite superior de integración (a) y presione Enter.'
	read( *, * ) b
	write( *, *) "Ingrese un entero para la variable function_type y presione Enter."
	read(*,*) function_type
	write( *, * ) '---------------------------------------'
			
	20 format (A15, F15.6)
	21 format (A5, I8, A10, F8.6)
	
	do k = 2, 8, 2
	
		n = 2**k
		m = n + 1
		h = ( b - a ) * ( 1._dp / n )
		exact_integ = -( f(b,function_type) - f(a,function_type) )
		m_gauss = 2**k
		
		
		write( *, 21 ) 'n = ', n, 'h  = ', h
		write( *, 20 ) 'I_exact       = ', exact_integ
		call trapez_integ ( m, a, b, h, trapez_num_integ, function_type )
		write( *, 20 ) 'I_trapezoidal = ', trapez_num_integ
		call errors_num_integrals(exact_integ, trapez_num_integ(1,1), rel_error, 2)
		write( *, 20 ) 'relative error = ', rel_error*100._dp
		
		call simpson_13_integ ( m, a, b, h, simps_13_num_integ, function_type )
		write( *, 20 ) 'I_simpson13   = ', simps_13_num_integ
		call errors_num_integrals(exact_integ, simps_13_num_integ(1,1), rel_error, 2)
		write( *, 20 ) 'relative error = ', rel_error*100._dp
		
		call simpson_38_integ ( m, a, b, h, simps_38_num_integ, function_type )
		write( *, 20 ) 'I_simpson38   = ', simps_38_num_integ
		call errors_num_integrals(exact_integ, simps_38_num_integ(1,1), rel_error, 2)
		write( *, 20 ) 'relative error = ', rel_error*100._dp
		
		call gauss_integ ( m_gauss, 0, a, b, gauss_num_integ, function_type)
		write( *, 20 ) 'I_gauss       = ', gauss_num_integ
		call errors_num_integrals(exact_integ, gauss_num_integ(1,1), rel_error, 2)
		write( *, 20 ) 'relative error = ', rel_error*100._dp
		
		write( *, * ) '---------------------------------------'
	end do
	
	return
end program num_integ_gauss


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
! gfortran -o num_integ_gauss ../../modules/module_presition.f90 ../../modules/module_functions.f90 ../../modules/module_gauss.f90 ../../modules/module_num_integrals.f90 ../../modules/module_numerical_error.f90 num_integ_gauss.f90 && ./num_integ_gauss
!----------------------------------------------------------
