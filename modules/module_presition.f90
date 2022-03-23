!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
! Modulo de precisión
!----------------------------------------------------------

module module_presition
	implicit none
	
	integer, parameter :: sp = selected_real_kind( p=6, r=37 )		! simple presicion (sp) class
	integer, parameter :: dp = selected_real_kind( p=15, r=307 )	! double presicion (dp) class
	
end module module_presition
