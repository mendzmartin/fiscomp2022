!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
!
!----------------------------------------------------------

program nombre

	implicit none

	! Explicit variables declaration
	integer, parameter :: sp = selected_real_kind( p=6, r=37 )		! simple presicion (sp) class
	integer, parameter :: dp = selected_real_kind( p=15, r=307 )	! double presicion (dp) class
	integer, parameter :: qp = selected_real_kind( p=33, r=4931 )	! quad presicion (dp) class
	
	real(sp) :: pi_acos1, pi_acos2, pi_acos3, pi_atan, pi_asin
	
	pi_atan = 4._sp*atan(1._sp)
	pi_acos1 = 2._sp*acos(0._sp)
	pi_acos2 = acos(-1._sp)
	pi_acos3 = 3._sp*acos(0.5_sp)
	pi_asin = 2._sp*asin(1._sp)
	
	!diff = abs(pi_atan - pi_acos)
	
	21 format (A15, F32.24)
	
	write(*,21) 'option 1 = ', pi_atan
	write(*,21) 'option 2 = ', pi_acos1
	write(*,21) 'option 3 = ', pi_acos2
	write(*,21) 'option 4 = ', pi_acos3
	write(*,21) 'option 5 = ', pi_asin
	!write(*,*) 'diff = ', diff
	
end program nombre


!----------------------------------------------------------
! REFERENCES
!----------------------------------------------------------
! https://stackoverflow.com/questions/2157920/why-define-pi-4atan1-d0
! http://computer-programming-forum.com/49-fortran/9088355fc7746aa5.htm
!
!----------------------------------------------------------

!----------------------------------------------------------
! COMPILATION RUL
!----------------------------------------------------------
! gfortran -g3 -Wall -Wextra -Wconversion nombre.f90 -o nombre && ./nombre
!
!
!----------------------------------------------------------
