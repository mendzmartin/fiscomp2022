!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
! Integración numérica.
! Puntos equiespaciados: Comparación Trapezoidal vs. Simpson
! 	programa sólo válido para la integral definida de la
!----------------------------------------------------------

program num_integ
	
	use presition_module 		! use module of presition
	use functions_module 		! use module of functions
	use num_integrals_module 	! use module of numerical integrals
	
	implicit none

	! Explicit variables declaration
	real(dp)					:: exact_integ				! valor de la integral exacta
	real(dp), dimension(1,1)	:: trapez_num_integ			! numerical integration with trapezoidal method
	real(dp), dimension(1,1) 	:: simps_13_num_integ		! numerical integration with simpson 1/3 method
	real(dp), dimension(1,1) 	:: simps_38_num_integ		! numerical integration with simpson 3/8 method
	integer(dp)					:: k						! potencia
	integer(dp) 				:: n, m 					! número de intervalos y número de puntos
	real(dp)					:: a,b						! límites de integración
	real(dp)					:: h						! paso de integración
	
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
		m = n + 1_dp 
		h = ( b - a ) * ( 1_dp / n )
		exact_integ = ( exp_func( b ) - exp_func(a) )
		
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
