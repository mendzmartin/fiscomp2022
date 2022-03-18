!----------------------------------------------------------
! Purpose:
!  Problem 2. Diferenciacion numérica evaluada en un valor
!  x particular utilizando la fórmula centrada en dos puntos
!----------------------------------------------------------

program diff_num

	implicit none

	! Explicit variables declaration
	integer, parameter :: dp = selected_real_kind( p=15, r=307 )	! double presition class
	real(dp) :: h, error, f_forw, f_back, diff_exact, diff_aprox, x	! reals of "dp" class
	integer(dp) :: i, i_start, i_end, delta_i						! integer of "dp" class
	integer(4) :: istat
	
	! user specific value
	write( *, * ) '########## ########## ########## ##########'
	write( *, * ) 'Hello, please input the point evaluation:'
	read( *, * ) x
	
	diff_exact = x*exp(x)	! exact differentiation

	open( 10, file = './result.dat', status = 'replace', iostat = istat )
	
	write(*,*) 'Input/Output file. istat = ', istat
	write( *, * ) '########## ########## ########## ##########'
	
	i_start	= 0
	i_end	= -2
	delta_i	= -1
	
	do i = 0, -20, -1	! var = start, end [,delta]
	
		h = 10d0**i												! real step = {1,0.1,0.01}
		f_forw = exp(x+h)										! forward function
		f_back = exp(x-h)										! backward function 
		diff_aprox = (f_forw - f_back)/(2*h)					! aprox differentiation
		error = abs((diff_exact - diff_aprox)/diff_exact)*100	! relative error
		write(10,*) -i, h, diff_exact, diff_aprox, error
	enddo
	
	close(10)
	
end program diff_num

! https://fortranwiki.org/fortran/show/Real+precision

!----------------------------------------------------------
! Regla de compilación
!----------------------------------------------------------
! gfortran -g3 -Wall -Wextra -Wconversion -o diff_num diff_num.f90 && rm result.dat && ./diff_num
