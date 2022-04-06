!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
! Modulo de precisión
!----------------------------------------------------------

module module_presition
	implicit none
	
	integer, parameter :: sp = selected_real_kind( p=6, r=37 )		! simple presicion (sp) class
	integer, parameter :: dp = selected_real_kind( p=15, r=307 )	! double presicion (dp) class
	integer, parameter :: qp = selected_real_kind( p=33, r=4931 )	! quad presicion (dp) class
	
end module module_presition

!----------------------------------------------------------
! REFERENCES
!----------------------------------------------------------
! https://fortranwiki.org/fortran/show/Real+precision
!----------------------------------------------------------
