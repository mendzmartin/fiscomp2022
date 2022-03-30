!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
!
!
!----------------------------------------------------------

program ODE_second_orden

	use module_presition 			! module for define presition
	use module_functions_1D 		! module for use 1D functions
	use module_functions_2D 		! module for use 2D functions
	use module_EDO_segundo_orden 	! module for resolve second order ODEs
	use module_numerical_error 		! module for estimate numerical errors

	implicit none
	
	! Explicit variables declaration
	integer(sp), parameter 	:: n = 8192					! points
	integer(sp) 			:: i,index_prev 			! index loop
	integer(sp) 			:: istat					! integer of simple presition
	real(dp)				:: a,b						! intervalo de validez para la solución aproximada
	real(dp)				:: h						! time step
	real(dp) 				:: y1_0, y2_0				! initial condition
	real(dp) 				:: rel_error_euler1 		! relative error Euler method
	real(dp) 				:: rel_error_RK2_hu1		! relative error 2nd order RK-Hung method
	real(dp) 				:: rel_error_RK2_mp1		! relative error 2nd order RK-Middle-Point method
	real(dp) 				:: rel_error_RK2_ra1		! relative error 2nd order RK-Raltson method
	real(dp) 				:: rel_error_RK4_cl1		! relative error 4th order RK-Classic method
	real(dp) 				:: rel_error_euler2 		! relative error Euler method
	real(dp) 				:: rel_error_RK2_hu2		! relative error 2nd order RK-Hung method
	real(dp) 				:: rel_error_RK2_mp2		! relative error 2nd order RK-Middle-Point method
	real(dp) 				:: rel_error_RK2_ra2		! relative error 2nd order RK-Raltson method
	real(dp) 				:: rel_error_RK4_cl2		! relative error 4th order RK-Classic method
	real(dp), dimension(n) 	:: y1_exact, y2_exact		! vector con solución exacta de la EDO
	real(dp), dimension(n) 	:: y1_euler, y2_euler 		! vector con solución aproximada de la EDO usando método de Euler
	real(dp), dimension(n) 	:: y1_RK2_hu, y2_RK2_hu		! aproximate solution vector of ODE using 2nd order RK-Hung method
	real(dp), dimension(n) 	:: y1_RK2_mp, y2_RK2_mp		! aproximate solution vector of ODE using 2nd order RK-Middle-Point method
	real(dp), dimension(n) 	:: y1_RK2_ra, y2_RK2_ra		! aproximate solution vector of ODE using 2nd order RK-Raltson method
	real(dp), dimension(n) 	:: y1_RK4_cl, y2_RK4_cl 	! aproximate soluction vector of ODE using 4th order RK-Classic method
	real(dp), dimension(n) 	:: x 						! variable vector
	real(dp), dimension(n)	:: engy_kin_exact, engy_pot_exact 	! exact kinetic and potential energies
	real(dp), dimension(n)	:: engy_kin_euler, engy_pot_euler 	! euler kinetic and potential energies
	real(dp), dimension(n)	:: engy_kin_RK2_hu, engy_pot_RK2_hu ! RK2_hu kinetic and potential energies
	real(dp), dimension(n)	:: engy_kin_RK2_mp, engy_pot_RK2_mp ! RK2_mp kinetic and potential energies
	real(dp), dimension(n)	:: engy_kin_RK2_ra, engy_pot_RK2_ra ! RK2_ra kinetic and potential energies
	real(dp), dimension(n)	:: engy_kin_RK4_cl, engy_pot_RK4_cl ! RK4_cl kinetic and potential energies
	
	! Datos pedidos al usuario
	write( *, * ) 'Ingrese el limite inferior de integración (a) y presione Enter.'
	read( *, * ) a
	write( *, * ) 'Ingrese el limite superior de integración (a) y presione Enter.'
	read( *, * ) b
	write( *, * ) 'Ingrese la condición inicial dy/dx para x=0 y presione Enter.'
	read( *, * ) y1_0
	write( *, * ) 'Ingrese la condición inicial (d^2y)/(dx^2) para x=0 y presione Enter.'
	read( *, * ) y2_0
	
	open( 10, file = './result_y1.dat', status = 'replace', iostat = istat )
	write(*,*) 'Input/Output file. istat = ', istat
	open( 11, file = './result_y2.dat', status = 'replace', iostat = istat )
	write(*,*) 'Input/Output file. istat = ', istat
	open( 12, file = './result_energies.dat', status = 'replace', iostat = istat )
	write(*,*) 'Input/Output file. istat = ', istat
	20 format (E12.4, E12.4, E12.4, E12.4, E12.4, E12.4, E12.4, E12.4, E12.4, E12.4, E12.4, E12.4)
	21 format (E12.4, E12.4, E12.4, E12.4, E12.4, E12.4)
	
	h = abs(b - a) * ( 1._dp / (n-1_sp) )
	
	x(1) = a
	y1_exact(1) = f_1D_HO(1._dp,1._dp,y1_0,y2_0,x(1), 1_sp) 	! f_1D_HO(m,k,y1_0,y2_0,x,function_type)
	y2_exact(1) = f_1D_HO(1._dp,1._dp,y1_0,y2_0,x(1), 2_sp) 	! f_1D_HO(m,k,y1_0,y2_0,x,function_type)
	
	engy_kin_exact = f_1D_HO(1._dp,1._dp,y1_0,y2_0,x(1), 4_sp)
	engy_pot_exact = f_1D_HO(1._dp,1._dp,y1_0,y2_0,x(1), 5_sp)
	
	call euler_method(n,a,b,y1_0,y2_0,y1_euler,y2_euler,0,1)	! calculamos solución con método de euler
	call RK2(1,n,a,b,y1_0,y2_0,y1_RK2_hu,y2_RK2_hu,0,1) 		! calculamos solución con método de RK2_hu
	call RK2(2,n,a,b,y1_0,y2_0,y1_RK2_mp,y2_RK2_mp,0,1) 		! calculamos solución con método de RK2_mp
	call RK2(3,n,a,b,y1_0,y2_0,y1_RK2_ra,y2_RK2_ra,0,1) 		! calculamos solución con método de RK2_ra
	call RK4(1,n,a,b,y1_0,y2_0,y1_RK4_cl,y2_RK4_cl,0,1)			! calculamos solución con método de RK4_mp
	
	call energies(n, 1._dp, 1._dp, y1_euler, y2_euler, engy_kin_euler, engy_pot_euler) 		! calculamos energias para euler
	call energies(n, 1._dp, 1._dp, y1_RK2_hu, y2_RK2_hu, engy_kin_RK2_hu, engy_pot_RK2_hu) 	! calculamos energias para RK2_hu
	call energies(n, 1._dp, 1._dp, y1_RK2_mp, y2_RK2_mp, engy_kin_RK2_mp, engy_pot_RK2_mp) 	! calculamos energias para RK2_mp
	call energies(n, 1._dp, 1._dp, y1_RK2_ra, y2_RK2_ra, engy_kin_RK2_ra, engy_pot_RK2_ra) 	! calculamos energias para RK2_ra
	call energies(n, 1._dp, 1._dp, y1_RK4_cl, y2_RK4_cl, engy_kin_RK4_cl, engy_pot_RK4_cl) 	! calculamos energias para RK4_cl
	
	call basic_errors_num(y1_exact(1), y1_euler(1), rel_error_euler1, 2) 	! calculamos error relativo método de euler
	call basic_errors_num(y1_exact(1), y1_RK2_hu(1), rel_error_RK2_hu1, 2) 	! calculamos error relativo método de RK2_hu
	call basic_errors_num(y1_exact(1), y1_RK2_mp(1), rel_error_RK2_mp1, 2) 	! calculamos error relativo método de RK2_mp
	call basic_errors_num(y1_exact(1), y1_RK2_ra(1), rel_error_RK2_ra1, 2) 	! calculamos error relativo método de RK2_ra
	call basic_errors_num(y1_exact(1), y1_RK4_cl(1), rel_error_RK4_cl1, 2) 	! calculamos error relativo método de RK4_cl
	
	call basic_errors_num(y2_exact(1), y2_euler(1), rel_error_euler2, 2) 	! calculamos error relativo método de euler
	call basic_errors_num(y2_exact(1), y2_RK2_hu(1), rel_error_RK2_hu2, 2) 	! calculamos error relativo método de RK2_hu
	call basic_errors_num(y2_exact(1), y2_RK2_mp(1), rel_error_RK2_mp2, 2) 	! calculamos error relativo método de RK2_mp
	call basic_errors_num(y2_exact(1), y2_RK2_ra(1), rel_error_RK2_ra2, 2) 	! calculamos error relativo método de RK2_ra
	call basic_errors_num(y2_exact(1), y2_RK4_cl(1), rel_error_RK4_cl2, 2) 	! calculamos error relativo método de RK4_cl
	
	
	write(10,20) 	x(1), y1_exact(1),&
					y1_euler(1),&
					y1_RK2_hu(1), y1_RK2_mp(1), y1_RK2_ra(1),&
					y1_RK4_cl(1),&
					rel_error_euler1,&
					rel_error_RK2_hu1, rel_error_RK2_mp1, rel_error_RK2_ra1,&
					rel_error_RK4_cl1
					
	write(11,20) 	x(1), y2_exact(1),&
					y2_euler(1),&
					y2_RK2_hu(1), y2_RK2_mp(1), y2_RK2_ra(1),&
					y2_RK4_cl(1),&
					rel_error_euler2,&
					rel_error_RK2_hu2, rel_error_RK2_mp2, rel_error_RK2_ra2,&
					rel_error_RK4_cl2
	
	write(12,21)	(engy_kin_exact(1)+engy_pot_exact(1)),&
					(engy_kin_euler(1)+engy_pot_euler(1)),&
					(engy_kin_RK2_hu(1)+engy_pot_RK2_hu(1)),&
					(engy_kin_RK2_mp(1)+engy_pot_RK2_mp(1)),&
					(engy_kin_RK2_ra(1)+engy_pot_RK2_ra(1)),&
					(engy_kin_RK4_cl(1)+engy_pot_RK4_cl(1))
	
	do i = 2, n, 1
		index_prev = i-1
		x(i) = x(index_prev) + h
		y1_exact(i) = f_1D_HO(1._dp,1._dp,y1_0,y2_0,x(i), 1_sp) 	! f_1D_HO(m,k,y1_0,y2_0,x,function_type)
		y2_exact(i) = f_1D_HO(1._dp,1._dp,y1_0,y2_0,x(i), 2_sp) 	! f_1D_HO(m,k,y1_0,y2_0,x,function_type)
		
		engy_kin_exact = f_1D_HO(1._dp,1._dp,y1_0,y2_0,x(i), 4_sp)
		engy_pot_exact = f_1D_HO(1._dp,1._dp,y1_0,y2_0,x(i), 5_sp)
		
		call basic_errors_num(y1_exact(i), y1_euler(i), rel_error_euler1, 2) 	! calculamos error relativo método de euler
		call basic_errors_num(y1_exact(i), y1_RK2_hu(i), rel_error_RK2_hu1, 2) 	! calculamos error relativo método de RK2
		call basic_errors_num(y1_exact(i), y1_RK2_mp(i), rel_error_RK2_mp1, 2) 	! calculamos error relativo método de RK2
		call basic_errors_num(y1_exact(i), y1_RK2_ra(i), rel_error_RK2_ra1, 2) 	! calculamos error relativo método de RK2
		call basic_errors_num(y1_exact(i), y1_RK4_cl(i), rel_error_RK4_cl1, 2) 	! calculamos error relativo método de RK4
		
		call basic_errors_num(y2_exact(i), y2_euler(i), rel_error_euler2, 2) 	! calculamos error relativo método de euler
		call basic_errors_num(y2_exact(i), y2_RK2_hu(i), rel_error_RK2_hu2, 2) 	! calculamos error relativo método de RK2
		call basic_errors_num(y2_exact(i), y2_RK2_mp(i), rel_error_RK2_mp2, 2) 	! calculamos error relativo método de RK2
		call basic_errors_num(y2_exact(i), y2_RK2_ra(i), rel_error_RK2_ra2, 2) 	! calculamos error relativo método de RK2
		call basic_errors_num(y2_exact(i), y2_RK4_cl(i), rel_error_RK4_cl2, 2) 	! calculamos error relativo método de RK4
		
		write(10,20) 	x(i), y1_exact(i),&
						y1_euler(i),&
						y1_RK2_hu(i), y1_RK2_mp(i), y1_RK2_ra(i),&
						y1_RK4_cl(i),&
						rel_error_euler1,&
						rel_error_RK2_hu1, rel_error_RK2_mp1, rel_error_RK2_ra1,&
						rel_error_RK4_cl1
						
		write(11,20) 	x(i), y2_exact(i),&
						y2_euler(i),&
						y2_RK2_hu(i), y2_RK2_mp(i), y2_RK2_ra(i),&
						y2_RK4_cl(i),&
						rel_error_euler2,&
						rel_error_RK2_hu2, rel_error_RK2_mp2, rel_error_RK2_ra2,&
						rel_error_RK4_cl2
						
		write(12,21)	(engy_kin_exact(i)+engy_pot_exact(i)),&
						(engy_kin_euler(i)+engy_pot_euler(i)),&
						(engy_kin_RK2_hu(i)+engy_pot_RK2_hu(i)),&
						(engy_kin_RK2_mp(i)+engy_pot_RK2_mp(i)),&
						(engy_kin_RK2_ra(i)+engy_pot_RK2_ra(i)),&
						(engy_kin_RK4_cl(i)+engy_pot_RK4_cl(i))
	end do
	
	close(10)
	close(11)
	close(12)
	
end program ODE_second_orden

subroutine energies(n, m, k, y1_vector, y2_vector, engy_kin, engy_pot)
	use module_presition 
	
	implicit none
	
	integer(sp), 			intent (in) :: n					! points 
	real(dp), 				intent(in) 	:: m, k 				! inertial mass & spring constant
	real(dp), dimension(n), intent(in) 	:: y1_vector, y2_vector ! position & velocity vectors
	real(dp), dimension(n), intent(out) 	:: engy_kin, engy_pot 	! kinetic and potential energies
	
	engy_kin = 0.5_dp*m*dot_product(y2_vector, y2_vector)
	engy_pot = 0.5_dp*k*dot_product(y1_vector, y1_vector)
	
end subroutine energies

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
! ./script_run.sh
!----------------------------------------------------------
