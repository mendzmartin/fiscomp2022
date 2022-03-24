!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
!
!
! Record of revisions
!  		Date 		Programmer
!		====		==========
!	20/03/2022		Martín Méndez
!----------------------------------------------------------

program nombre

	implicit none

	! Explicit variables declaration
	integer, parameter :: sp = selected_real_kind( p=6, r=37 )		! simple presicion (sp) class
	integer, parameter :: dp = selected_real_kind( p=15, r=307 )	! double presicion (dp) class
	
	call subroutine_name( argument_list )
	
end program nombre


!----------------------------------------------------------
! PURPOSE
!
!
!----------------------------------------------------------
subroutine subroutine_name ( argument_list )
	! Declaration section
	
	! Execution section
	
	return
end subroutine simpson_integ


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
! gfortran -g3 -Wall -Wextra -Wconversion nombre.f90 -o nombre && ./nombre
!
!
!----------------------------------------------------------
