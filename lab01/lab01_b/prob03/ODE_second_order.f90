!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
! Algorito para resolver la EDO de 2do orden del oscilador
!  armónico libre (sin fuerza externa), utilizando los mé-
!  todos de Euler, Runge Kutta de 2do orden y Runge Kutta
!  de 4to orden.
!----------------------------------------------------------

program ODE_second_orden_v2

	use module_presition 			! module for define presition
	use module_functions_1D 		! module for use 1D functions
	use module_functions_2D 		! module for use 2D functions
	use module_EDO_segundo_orden 	! module for resolve second order ODEs
	use module_numerical_error 		! module for estimate numerical errors

	implicit none
	
	! Explicit variables declaration
	integer(sp)							:: k						! pow
	integer(sp)							:: n						! points
	integer(sp) 						:: i,index_prev 			! index loop
	integer(sp) 						:: istat					! integer of simple presition
	real(dp)							:: a,b						! intervalo de validez para la solución aproximada
	real(dp)							:: h						! time step
	real(dp) 							:: y1_0, y2_0				! initial condition
	real(dp) 							:: rel_error_euler1 		! relative error Euler method
	real(dp) 							:: rel_error_RK2_hu1		! relative error 2nd order RK-Hung method
	real(dp) 							:: rel_error_RK2_mp1		! relative error 2nd order RK-Middle-Point method
	real(dp) 							:: rel_error_RK2_ra1		! relative error 2nd order RK-Raltson method
	real(dp) 							:: rel_error_RK4_cl1		! relative error 4th order RK-Classic method
	real(dp) 							:: rel_error_euler2 		! relative error Euler method
	real(dp) 							:: rel_error_RK2_hu2		! relative error 2nd order RK-Hung method
	real(dp) 							:: rel_error_RK2_mp2		! relative error 2nd order RK-Middle-Point method
	real(dp) 							:: rel_error_RK2_ra2		! relative error 2nd order RK-Raltson method
	real(dp) 							:: rel_error_RK4_cl2		! relative error 4th order RK-Classic method
	real(dp), dimension(:), allocatable :: y1_exact, y2_exact		! vector con solución exacta de la EDO
	real(dp), dimension(:), allocatable :: y1_euler, y2_euler 		! vector con solución aproximada de la EDO usando método de Euler
	real(dp), dimension(:), allocatable :: y1_RK2_hu, y2_RK2_hu		! aproximate solution vector of ODE using 2nd order RK-Hung method
	real(dp), dimension(:), allocatable :: y1_RK2_mp, y2_RK2_mp		! aproximate solution vector of ODE using 2nd order RK-Middle-Point method
	real(dp), dimension(:), allocatable :: y1_RK2_ra, y2_RK2_ra		! aproximate solution vector of ODE using 2nd order RK-Raltson method
	real(dp), dimension(:), allocatable :: y1_RK4_cl, y2_RK4_cl 	! aproximate soluction vector of ODE using 4th order RK-Classic method
	real(dp), dimension(:), allocatable :: x 						! variable vector
	real(dp), dimension(:), allocatable	:: engy_kin_exact, engy_pot_exact 	! exact kinetic and potential energies
	real(dp), dimension(:), allocatable	:: engy_kin_euler, engy_pot_euler 	! euler kinetic and potential energies
	real(dp), dimension(:), allocatable	:: engy_kin_RK2_hu, engy_pot_RK2_hu ! RK2_hu kinetic and potential energies
	real(dp), dimension(:), allocatable	:: engy_kin_RK2_mp, engy_pot_RK2_mp ! RK2_mp kinetic and potential energies
	real(dp), dimension(:), allocatable	:: engy_kin_RK2_ra, engy_pot_RK2_ra ! RK2_ra kinetic and potential energies
	real(dp), dimension(:), allocatable	:: engy_kin_RK4_cl, engy_pot_RK4_cl ! RK4_cl kinetic and potential energies
	
	! Datos pedidos al usuario
	write( *, * ) 'Ingrese el limite inferior de integración (a) y presione Enter.'
	read( *, * ) a
	write( *, * ) 'Ingrese el limite superior de integración (a) y presione Enter.'
	read( *, * ) b
	write( *, * ) 'Ingrese la condición inicial dy/dx para x=0 y presione Enter.'
	read( *, * ) y1_0
	write( *, * ) 'Ingrese la condición inicial (d^2y)/(dx^2) para x=0 y presione Enter.'
	read( *, * ) y2_0
	
	open( 10, file = './result_y1.dat', status = 'replace', action = 'write', iostat = istat )
	write(*,*) 'Input/Output file. istat = ', istat
	open( 11, file = './result_y2.dat', status = 'replace', action = 'write', iostat = istat )
	write(*,*) 'Input/Output file. istat = ', istat
	open( 12, file = './result_energies.dat', status = 'replace', action = 'write', iostat = istat )
	write(*,*) 'Input/Output file. istat = ', istat
	20 format (I10, E12.4, E12.4, E12.4, E12.4, E12.4, E12.4, E12.4, E12.4, E12.4, E12.4, E12.4, E12.4)
	21 format (I10, E12.4, E12.4, E12.4, E12.4, E12.4, E12.4, E12.4)
	
	if (istat == 0_sp) then
		do k = 1, 17, 1 ! (to obtain h [1, 10^-6])
		
			n = 2_dp**k ! n = {2,4,8,16,32,64,128,256,512,1024,...}
			
			h = abs(b - a) * ( 1._dp / (n-1_sp) )
			
			allocate(x(n))
			x(1) = a
			allocate(y1_exact(n),y2_exact(n))
			y1_exact(1) = f_1D_HO(1._dp,1._dp,y1_0,y2_0,x(1), 1_sp) 	! f_1D_HO(m,k,y1_0,y2_0,x,function_type)
			y2_exact(1) = f_1D_HO(1._dp,1._dp,y1_0,y2_0,x(1), 2_sp) 	! f_1D_HO(m,k,y1_0,y2_0,x,function_type)
			
			allocate(engy_kin_exact(n),engy_pot_exact(n))
			engy_kin_exact = f_1D_HO(1._dp,1._dp,y1_0,y2_0,x(1), 4_sp)
			engy_pot_exact = f_1D_HO(1._dp,1._dp,y1_0,y2_0,x(1), 5_sp)
			
			allocate(y1_euler(n),y2_euler(n))
			call euler_method(n,a,b,y2_0,y1_0,y2_euler,y1_euler,1,1)	! calculamos solución con método de euler
			allocate(y1_RK2_hu(n),y2_RK2_hu(n))
			call RK2(1,n,a,b,y2_0,y1_0,y2_RK2_hu,y1_RK2_hu,1,1) 		! calculamos solución con método de RK2_hu
			allocate(y1_RK2_mp(n),y2_RK2_mp(n))
			call RK2(2,n,a,b,y2_0,y1_0,y2_RK2_mp,y1_RK2_mp,1,1) 		! calculamos solución con método de RK2_mp
			allocate(y1_RK2_ra(n),y2_RK2_ra(n))
			call RK2(3,n,a,b,y2_0,y1_0,y2_RK2_ra,y1_RK2_ra,1,1) 		! calculamos solución con método de RK2_ra
			allocate(y1_RK4_cl(n),y2_RK4_cl(n))
			call RK4(1,n,a,b,y2_0,y1_0,y2_RK4_cl,y1_RK4_cl,1,1)			! calculamos solución con método de RK4_mp
			
			! calculamos energias para euler
			allocate(engy_kin_euler(n),engy_pot_euler(n))
			call energies(1._dp, 1._dp, y1_euler(1), y2_euler(1), engy_kin_euler(1), engy_pot_euler(1))
			! calculamos energias para RK2_hu
			allocate(engy_kin_RK2_hu(n),engy_pot_RK2_hu(n))
			call energies(1._dp, 1._dp, y1_RK2_hu(1), y2_RK2_hu(1), engy_kin_RK2_hu(1), engy_pot_RK2_hu(1))
			! calculamos energias para RK2_mp
			allocate(engy_kin_RK2_mp(n),engy_pot_RK2_mp(n))
			call energies(1._dp, 1._dp, y1_RK2_mp(1), y2_RK2_mp(1), engy_kin_RK2_mp(1), engy_pot_RK2_mp(1))
			! calculamos energias para RK2_ra
			allocate(engy_kin_RK2_ra(n),engy_pot_RK2_ra(n))
			call energies(1._dp, 1._dp, y1_RK2_ra(1), y2_RK2_ra(1), engy_kin_RK2_ra(1), engy_pot_RK2_ra(1))
			! calculamos energias para RK4_cl
			allocate(engy_kin_RK4_cl(n),engy_pot_RK4_cl(n))
			call energies(1._dp, 1._dp, y1_RK4_cl(1), y2_RK4_cl(1), engy_kin_RK4_cl(1), engy_pot_RK4_cl(1))
			
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
			
			
			write(10,20) 	n, x(1), y1_exact(1),&
							y1_euler(1),&
							y1_RK2_hu(1), y1_RK2_mp(1), y1_RK2_ra(1),&
							y1_RK4_cl(1),&
							rel_error_euler1,&
							rel_error_RK2_hu1, rel_error_RK2_mp1, rel_error_RK2_ra1,&
							rel_error_RK4_cl1
							
			write(11,20) 	n, x(1), y2_exact(1),&
							y2_euler(1),&
							y2_RK2_hu(1), y2_RK2_mp(1), y2_RK2_ra(1),&
							y2_RK4_cl(1),&
							rel_error_euler2,&
							rel_error_RK2_hu2, rel_error_RK2_mp2, rel_error_RK2_ra2,&
							rel_error_RK4_cl2
			
			write(12,21)	n, x(1), (engy_kin_exact(1)+engy_pot_exact(1)),&
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
				
				! calculamos energias para euler
				call energies(1._dp, 1._dp, y1_euler(i), y2_euler(i), engy_kin_euler(i), engy_pot_euler(i))
				! calculamos energias para RK2_hu
				call energies(1._dp, 1._dp, y1_RK2_hu(i), y2_RK2_hu(i), engy_kin_RK2_hu(i), engy_pot_RK2_hu(i))
				! calculamos energias para RK2_mp
				call energies(1._dp, 1._dp, y1_RK2_mp(i), y2_RK2_mp(i), engy_kin_RK2_mp(i), engy_pot_RK2_mp(i))
				! calculamos energias para RK2_ra
				call energies(1._dp, 1._dp, y1_RK2_ra(i), y2_RK2_ra(i), engy_kin_RK2_ra(i), engy_pot_RK2_ra(i))
				! calculamos energias para RK4_cl
				call energies(1._dp, 1._dp, y1_RK4_cl(i), y2_RK4_cl(i), engy_kin_RK4_cl(i), engy_pot_RK4_cl(i))
				
				write(10,20) 	n, x(i), y1_exact(i),&
								y1_euler(i),&
								y1_RK2_hu(i), y1_RK2_mp(i), y1_RK2_ra(i),&
								y1_RK4_cl(i),&
								rel_error_euler1,&
								rel_error_RK2_hu1, rel_error_RK2_mp1, rel_error_RK2_ra1,&
								rel_error_RK4_cl1
								
				write(11,20) 	n, x(i), y2_exact(i),&
								y2_euler(i),&
								y2_RK2_hu(i), y2_RK2_mp(i), y2_RK2_ra(i),&
								y2_RK4_cl(i),&
								rel_error_euler2,&
								rel_error_RK2_hu2, rel_error_RK2_mp2, rel_error_RK2_ra2,&
								rel_error_RK4_cl2
								
				write(12,21)	n, x(i), (engy_kin_exact(i)+engy_pot_exact(i)),&
								(engy_kin_euler(i)+engy_pot_euler(i)),&
								(engy_kin_RK2_hu(i)+engy_pot_RK2_hu(i)),&
								(engy_kin_RK2_mp(i)+engy_pot_RK2_mp(i)),&
								(engy_kin_RK2_ra(i)+engy_pot_RK2_ra(i)),&
								(engy_kin_RK4_cl(i)+engy_pot_RK4_cl(i))
			end do
			
			deallocate(y1_exact,y2_exact,x)
			deallocate(y1_euler,y2_euler)
			deallocate(y1_RK2_hu,y2_RK2_hu)
			deallocate(y1_RK2_mp,y2_RK2_mp)
			deallocate(y1_RK2_ra,y2_RK2_ra)
			deallocate(y1_RK4_cl,y2_RK4_cl)
			
			deallocate(engy_kin_exact,engy_pot_exact)
			deallocate(engy_kin_euler,engy_pot_euler)
			deallocate(engy_kin_RK2_hu,engy_pot_RK2_hu)
			deallocate(engy_kin_RK2_mp,engy_pot_RK2_mp)
			deallocate(engy_kin_RK2_ra,engy_pot_RK2_ra)
			deallocate(engy_kin_RK4_cl,engy_pot_RK4_cl)
			
		end do
		
		close(10)
		close(11)
		close(12)
	else
		write(*,*) 'Error opening file. iostat != 0'
	end if
end program ODE_second_orden_v2

!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
! Subrutina que sirve para calcular las energías cinética y
!  potencial del oscilador armónico unidimensional.
!----------------------------------------------------------
subroutine energies(m, k, y1, y2, engy_kin, engy_pot)
	use module_presition 
	
	implicit none
	
	real(dp), 		intent(in) 	:: m, k 				! inertial mass & spring constant
	real(dp), 		intent(in) 	:: y1, y2 				! position & velocity
	real(dp), 		intent(out) :: engy_kin, engy_pot 	! kinetic and potential energies
	
	engy_kin = 0.5_dp*m*(y2**2)
	engy_pot = 0.5_dp*k*(y1**2)
	
end subroutine energies

!----------------------------------------------------------
! REFERENCES
!----------------------------------------------------------
!
!
!----------------------------------------------------------

!----------------------------------------------------------
! COMPILATION RUL
!----------------------------------------------------------
! ./script_run.sh
!----------------------------------------------------------
