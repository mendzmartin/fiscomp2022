!----------------------------------------------------------
! Las siguientes expresiones son legales o ilegales en
!  Fortran 90? Si son legales cuál es su resultado? Si son
!  ilegales, qué hay de malo en ellas?
!----------------------------------------------------------
program prob01
	
	use module_presition
	
    implicit none
    
    real(dp) 	:: res 		! variable to save de aritmetical result in double presition
    integer(sp) :: istat 	! integer of simple presition
    
    open( unit = 10, file = 'results.dat', status = 'replace', iostat = istat )
	write(*,*) 'Input/Output file. istat = ', istat
	20 format (A15, F12.4)
    
    res = 37_dp/3_dp
    write(10,20) '37/3 = ', res
    
	res = 37_dp + 17_dp/3_dp
    write(10,20) '37+17/3 = ', res
    
	res = 28_dp/3_dp/4_dp
    write(10,20) '28/3/4 = ', res
	
	res = (28_dp/3_dp)/4_dp
    write(10,20) '(28/3)/4 = ', res
    
!	res = 28_dp/(3_dp/4_dp)
!	write(10,20) '28/(3/4) = ', res
    write(10,*) '28/(3/4)  =  Error: Division by zero'
    
	res = -3._dp**4._dp/2._dp
    write(10,20) '-3.**4./2. = ', res
    
	res = 3._dp**(-4._dp/2._dp)
    write(10,20) '3.**(-4./2.) = ', res
    
	res = 4._dp**-3_dp
    write(10,20) '4.**-3 = ', res
    
	res = 2._dp**2._dp**3._dp
    write(10,20) '2.**2.**3. = ', res
    
	res = 2._dp**(-2._dp)
    write(10,20) '2.**(-2.) = ', res
    
	res = (-2_dp)**2_dp
    write(10,20) '(-2)**2 = ', res

!	res = (-2._dp)**(-2.2_dp)
!	write(10,20) '(-2.)**(-2.2) = ', res
    write(10,*) '(-2.)**(-2.2)  =  Error: Raising a negative REAL at expresion to a REAL power is prohibited'
    
    
	res = (-2._dp)**nint(-2.2_dp)
    write(10,20) '(-2.)**nint(-2.2) = ', res
    
	res = 1_dp+1_dp/4_dp
    write(10,20) '1+1/4 = ', res
    
	res = 1._dp+1_dp/4_dp
    write(10,20) '1.+1/4 = ', res
    
	res = 1_dp*1._dp/4_dp
    write(10,20) '1*1./4 = ', res
    
    close(10)
    
    return
end program prob01
!----------------------------------------------------------
! Regla de compilación
!----------------------------------------------------------
! gfortran -g3 -Wall -Wextra -Wconversion -o probXX prob01.f90 && ./probXX

! gfortran ../../modules/module_presition.f90 -o prob01 prob01.f90 && ./prob01
