!----------------------------------------------------------
! Las siguientes expresiones son legales o ilegales en
!  Fortran 90? Si son legales cuál es su resultado? Si son
!  ilegales, qué hay de malo en ellas?
!----------------------------------------------------------
program prob01
    implicit none
    
    integer, parameter :: dp = kind( 1.0d0 ) ! double presition
    real(dp) :: res
    
    !res = 37/3						! prob01_01.o
	!res = 37 + 17 / 3				! prob01_02.o
	!res = 28 / 3 / 4				! prob01_03.o
	res = ( 28 / 3 ) / 4			! prob01_04.o
	!res = 28 / (3 / 4)				! prob01_05.o
	!res = -3. ** 4. / 2.			! prob01_06.o
	!res = -3. ** 4. / 2.			! prob01_07.o
	!res = 4. ** -3					! prob01_08.o
	!res = 2. ** 2. ** 3.			! prob01_09.o
	!res = 2. ** (-2.)				! prob01_10.o
	!res = (-2) ** 2				! prob01_11.o
	!res = (-2.) ** (-2.2)			! prob01_12.o
	!res = (-2.) ** NINT (-2.2)		! prob01_13.o
	!res = 1 + 1/4					! prob01_14.o
	!res = 1. + 1/4					! prob01_15.o
	!res = 1 * 1./4					! prob01_16.o
	
    print *, "result = ", res
    
    return
end program prob01
!----------------------------------------------------------
! Regla de compilación
!----------------------------------------------------------
! gfortran -g3 -Wall -Wextra -Wconversion -o probXX prob01.f90 && ./probXX
