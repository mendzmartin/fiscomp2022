program double_pendulum
	use module_presition
	use module_EDO_segundo_orden
    use module_double_pendulum
	use module_fourier_transform

	implicit none

	integer(sp), 	parameter 		:: n=500 !1048576								! points number
	real(dp), 		dimension(n) 	:: y1_RK4, y2_RK4, y3_RK4, y4_RK4 			! solutions from RK4 method
    real(dp),       dimension(n)    :: t                               			 ! time vector
	real(dp), 		parameter 		:: ti=0._dp, tf=450._dp 					! start and end times
	real(dp) 						:: y1_0,y2_0,y3_0,y4_0 						! initial conditions
	real(dp) 						:: h 										! step
	integer(sp) 					:: i,index_prev,k							! loop variables
	integer(sp) 					:: istat 									! control variable
    real(dp)                        :: a,b,c                            		! parameters 
    real(dp)                        :: kin_energy,pot_energy,tot_energy,rest 	! kinetic, potential and total energies
    
	real(dp) 						:: h_ic
	real(dp) 						:: tot_energy_exact, epsilon_energy
	real(dp) 						:: y3_0_min, y3_0_max
    integer(sp), parameter 			:: num_ic=100000
	logical 						:: control

	integer(sp), parameter 					:: n_st = 250 !524288 							! points state of steady state
	integer(sp) 							:: n_new
	real(dp), dimension(:), allocatable 	:: y1_RK4_st, y2_RK4_st 					! RK4 solutions of steady state
	real(dp), dimension(:), allocatable		:: ft_real_1, ft_aimag_1, sde_1, omega_1 	! fourier transform results
	real(dp), dimension(:), allocatable		:: ft_real_2, ft_aimag_2, sde_2, omega_2 	! fourier transform results


    a = (1._dp/3._dp) 	! alpha = m2/m1
	b = 0.5_dp		 	! beta  = l2/l1
	c = 0.5_dp			! gamma = g/l1

	y1_0 = 0._dp
	y2_0 = 0._dp
	y4_0 = 0._dp
	
	! Calculamos la condición inicial y3_0 para obtener la energía buscada
	! Como RK4 conseva la energía podemos calcularla en cualquier intervalo
	! 	temporal, en particular la calculamos para la condición inicial.

	y3_0_min = 0.0865_dp
	y3_0_max = 0.0867_dp
	
	h_ic = abs(y3_0_max-y3_0_min)/(real(num_ic,dp)-1._dp)
	
	tot_energy_exact = -0.745_dp
	epsilon_energy = 1e-10_dp

	control = .false.

	do i = 0, (num_ic-1)
		y3_0 = y3_0_min + h_ic*i
		call energy_dble_pendulum (y1_0, y2_0, y3_0,y4_0, a, b, c, kin_energy, pot_energy)
		tot_energy = (kin_energy + pot_energy)
		rest = abs(tot_energy - tot_energy_exact)
		if (rest <= epsilon_energy) then
			write(*,*) i, y3_0
			control = .true.
			exit
		end if
	end do

	if (control .eqv. .false.) then
		write(*,*) 'No se cumple la condicion para ningun valor'
		write(*,*) 'total energy = ', tot_energy
	else
		write(*,*) 'Se cumple la condicion para el valor y3_0 = ', y3_0
		write(*,*) 'total energy = ', tot_energy
	end if

	!subroutine RK4_four_eq(RK4_type, n, a, b, y1_0, y2_0, y3_0,y4_0, y1_RK4, y2_RK4, y3_RK4, y4_RK4, function_type, input_type)
	call RK4_four_eq(1_sp,n,ti,tf,y1_0,y2_0,y3_0,y4_0,y1_RK4,y2_RK4,y3_RK4,y4_RK4,1_sp,1_sp)

	n_new = (n-n_st)+1_sp
	allocate(y1_RK4_st(n_new), y2_RK4_st(n_new))
	! descartamos regimen )ransitorio de l)s soluciones
	y1_RK4_st = y1_RK4(n_st:n)
	y2_RK4_st = y2_RK4(n_st:n)


    t(1) = ti

    h = abs( tf - ti ) * ( 1._dp / ( real(n,dp) - 1._dp ) )
	do i= 2, n
        index_prev = i-1
		t(i) = t(1) + h*index_prev
	end do

	! REALIZAMOS LAS TRANSFORMADAS DE FOURIER DE LAS DOS SOLUCIONES (DESCARTANDO EL TRANSITORIO)

	allocate(ft_real_1(n_new+1_sp), ft_aimag_1(n_new+1_sp), sde_1(n_new+1_sp), omega_1(n_new+1_sp))
	allocate(ft_real_2(n_new+1_sp), ft_aimag_2(n_new+1_sp), sde_2(n_new+1_sp), omega_2(n_new+1_sp))
	! subroutine fft_forware_r2c(N, t_start, t_step, func, ft_real, ft_aimag, sde, omega)
	call fft_forware_r2c(n_new, t(n_st), t(n), y1_RK4_st, ft_real_1, ft_aimag_1, sde_1, omega_1)
	call fft_forware_r2c(n_new, t(n_st), t(n), y2_RK4_st, ft_real_2, ft_aimag_2, sde_2, omega_2)

	20 format (E14.6, E14.6, E14.6, E14.6)

	open( 10, file = './result_ft_01.dat', status = 'replace', action = 'write', iostat = istat )
	write(*,*) 'istat of result_ft_01.dat = ', istat
	open( 11, file = './result_ft_02.dat', status = 'replace', action = 'write', iostat = istat )
	write(*,*) 'istat of result_ft_02.dat = ', istat

	do k = -n_new/2, n_new/2
		if ( k < 0 ) then
			write(10,20) omega_1(-k), ft_real_1(-k), ft_aimag_1(-k), sde_1(-k)
			write(11,20) omega_2(-k), ft_real_2(-k), ft_aimag_2(-k), sde_2(-k)
		else
			write(10,20) omega_1(n_new/2+1+k), ft_real_1(n_new/2+1+k), ft_aimag_1(n_new/2+1+k), sde_1(n_new/2+1+k)
			write(11,20) omega_2(n_new/2+1+k), ft_real_2(n_new/2+1+k), ft_aimag_2(n_new/2+1+k), sde_2(n_new/2+1+k)
		end if
	end do

	! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	y3_0 = sqrt(1.125_dp)

	!subroutine RK4_four_eq(RK4_type, n, a, b, y1_0, y2_0, y3_0,y4_0, y1_RK4, y2_RK4, y3_RK4, y4_RK4, function_type, input_type)
	call RK4_four_eq(1_sp,n,ti,tf,y1_0,y2_0,y3_0,y4_0,y1_RK4,y2_RK4,y3_RK4,y4_RK4,1_sp,1_sp)

	! descartamos regimen )ransitorio de l)s soluciones
	y1_RK4_st = y1_RK4(n_st:n)
	y2_RK4_st = y2_RK4(n_st:n)

	! REALIZAMOS LAS TRANSFORMADAS DE FOURIER DE LAS DOS SOLUCIONES (DESCARTANDO EL TRANSITORIO)

	! subroutine fft_forware_r2c(N, t_start, t_step, func, ft_real, ft_aimag, sde, omega)
	call fft_forware_r2c(n_new, t(n_st), t(n), y1_RK4_st, ft_real_1, ft_aimag_1, sde_1, omega_1)
	call fft_forware_r2c(n_new, t(n_st), t(n), y2_RK4_st, ft_real_2, ft_aimag_2, sde_2, omega_2)

	open( 12, file = './result_ft_01_E0.dat', status = 'replace', action = 'write', iostat = istat )
	write(*,*) 'istat of result_ft_01_E0.dat = ', istat
	open( 13, file = './result_ft_02_E0.dat', status = 'replace', action = 'write', iostat = istat )
	write(*,*) 'istat of result_ft_02_E0.dat = ', istat

	do k = -n_new/2, n_new/2
		if ( k < 0 ) then
			write(12,20) omega_1(-k), ft_real_1(-k), ft_aimag_1(-k), sde_1(-k)
			write(13,20) omega_2(-k), ft_real_2(-k), ft_aimag_2(-k), sde_2(-k)
		else
			write(12,20) omega_1(n_new/2+1+k), ft_real_1(n_new/2+1+k), ft_aimag_1(n_new/2+1+k), sde_1(n_new/2+1+k)
			write(13,20) omega_2(n_new/2+1+k), ft_real_2(n_new/2+1+k), ft_aimag_2(n_new/2+1+k), sde_2(n_new/2+1+k)
		end if
	end do

	deallocate(y1_RK4_st, y2_RK4_st)
	deallocate(ft_real_1, ft_aimag_1, sde_1, omega_1)
	deallocate(ft_real_2, ft_aimag_2, sde_2, omega_2)

	close(10)
	close(11)
	close(12)
	close(13)


	! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end program double_pendulum

! OJO! NOTAR QUE LOS RESULTADOS VARIAN MUCHO RESPECTO A LA ELECCIÓN DEL NÚMERO "n"
! La variable tot_energy está

! https://web.mit.edu/jorloff/www/chaosTalk/double-pendulum/double-pendulum-en.html

