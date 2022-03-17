!----------------------------------------------------------
! Purpose:
!  Problem 2. Diferenciacion numérica
!----------------------------------------------------------

program diff_num

	implicit none

	integer, parameter :: dp = selected_real_kind( p=15, r=307 )	! double presition class
	real(dp) :: h, error						! real de clase dp
	real(dp) :: f, f_forw, f_back
	real(dp) :: diff_exact, diff_aprox
	real(dp) :: x
	integer(dp) :: i
	
	x = 1 ! valor
	
	do i = 0, 2, 1 ! var = start, stop [,step]
	
		h = 10**i
		f = exp(x) ! funcion
		f_forw = exp(x+h) ! funcion forward
		f_back = exp(x-h) ! funcino backward 
		diff_exact = x*f ! diferenciacion exacta
		diff_aprox = (f_forw - f_back)/(2*h) ! diferecnciacion aproximada
		error = abs((diff_exact - diff_aprox)/diff_exact)*100
		write(*,*) error
	enddo
	
	return
end program diff_num

!----------------------------------------------------------
! Regla de compilación
!----------------------------------------------------------
! gfortran -g3 -Wall -Wextra -Wconversion -o probXX prob01.f90 && ./probXX
