program hola

use module_presition
use module_functions_2D
use module_EDO_segundo_orden

implicit none

integer(sp), parameter :: n=1000

real(dp), dimension(n)	:: y1_RK4, y2_RK4, y3_RK4, y4_RK4
real(dp), parameter :: ti=0._dp, tf=450._dp 
real(dp) :: y1_0,y2_0,y3_0,y4_0,h
integer(sp) :: i, istat

y1_0=0._dp
y2_0=0._dp
y3_0=sqrt(1.125_dp)
y4_0=0._dp

h = abs(tf - ti) * ( 1._dp / (n-1_sp) )

!call RK4_four_eq(RK4_type, n, a, b, y1_0, y2_0, y3_0,y4_0, y1_RK4, y2_RK4, y3_RK4, y4_RK4, function_type, input_type)
call RK4_four_eq(1_sp, n, ti, tf, y1_0, y2_0, y3_0,y4_0, y1_RK4, y2_RK4, y3_RK4, y4_RK4, 1_sp, 1_sp)

open( 10, file = './result.dat', status = 'replace', action = 'write', iostat = istat )
write(*,*) 'Input/Output file. istat = ', istat

do i= 1, n
	write(10,*) real((i-1),dp)*h, y1_RK4(i),y2_RK4(i),y3_RK4(i),y4_RK4(i)
end do

close(10)

end program hola
