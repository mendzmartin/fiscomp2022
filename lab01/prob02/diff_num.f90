!----------------------------------------------------------
! PURPOSE
!  Problem 2. Diferenciacion numérica evaluada en un valor
!  x particular utilizando la fórmula centrada en dos puntos
!----------------------------------------------------------

program diff_num

	implicit none

	! Explicit variables declaration
	integer, parameter :: sp = selected_real_kind( p=6, r=37 )		! simple presicion (sp) class
	integer, parameter :: dp = selected_real_kind( p=15, r=307 )	! double presicion (dp) class
	real(dp) :: h, h_opt, x0, error, epsilon_machine				! reals of "dp" class
	real(dp) :: first_diff_exact, first_diff_aprox
	real(dp) :: f_forw, f_back
	integer(sp) :: i, i_start, i_end, delta_i						! integer of "sp" class
	integer(sp) :: istat											! integer of simple presition
	
	! user specific value
	write( *, * ) '########## ########## ########## ##########'
	write( *, * ) 'Hello, please input the point evaluation:'
	read( *, * ) x0
	
	first_diff_exact = exp(x0)			! first order exact differentiation
!	third_diff_exact = first_diff_exact	! third order exact differentiation

	!CALCULO DEL PASO OPTIMO
	
	! epsilon(x) returns the smallest number E of the same kind
	!  as x such that 1 + E > 1.
	epsilon_machine = epsilon(x0)
	
	h_opt = (24d0*epsilon_machine*(1/first_diff_exact))**(1d0/3d0)
	write(*,*) 'paso optimo = ', h_opt
	!30 format (E15.2)
	
	open( 10, file = './result.dat', status = 'new', iostat = istat )
	write(*,*) 'Input/Output file. istat = ', istat
	20 format (E15.2, E15.6)
	
	!CALCULO DE LA DIFERENCIACION NUMÉRICA
	!CALCULO DEL ERROR RELATIVO
	
	i_start	= 0
	i_end	= -10
	delta_i	= -1
	
	do i = i_start, i_end, delta_i	! var = start, end [,delta]
	
		h = 10d0**i																! real step = {1,0.1,0.01}
		f_forw = exp(x0+h)														! forward function
		f_back = exp(x0-h)														! backward function 
		first_diff_aprox = (f_forw - f_back)*(0.5d0/h)							! aprox differentiation
		error = abs((first_diff_exact - first_diff_aprox)*(1/first_diff_exact))	! relative error
!		error_porcen = error*100	! relative error
!		write(10,*) -i, h, diff_exact, diff_aprox, error
		write(10,20) h, error
	enddo
	
	f_forw = exp(x0+h_opt)														! forward function
	f_back = exp(x0-h_opt)														! backward function 
	first_diff_aprox = (f_forw - f_back)*(0.5d0/h_opt)							! aprox differentiation
	error = abs((first_diff_exact - first_diff_aprox)*(1/first_diff_exact))
	write(10,20) h_opt, error
	
	close(10)
	

	
end program diff_num

!----------------------------------------------------------
! REFERENCES
!----------------------------------------------------------
! https://fortranwiki.org/fortran/show/Real+precision
! https://gcc.gnu.org/onlinedocs/gfortran/EPSILON.html
! Landau pp. 113 (assessment: error analysis)

!----------------------------------------------------------
! COMPILATION RUL
!----------------------------------------------------------
! gfortran -g3 -Wall -Wextra -Wconversion -o diff_num diff_num.f90 && rm result.dat && ./diff_num && cat result.dat
! gfortran -o diff_num diff_num.f90 && rm result.dat && ./diff_num
