program double_pendulum
	use module_presition
	use module_EDO_segundo_orden
    use module_double_pendulum

	implicit none

	integer(sp), 	parameter 		:: n=1048576						! points number
	real(dp), 		dimension(n) 	:: y1_RK4, y2_RK4, y3_RK4, y4_RK4 	! solutions from RK4 method
    real(dp),       dimension(n)    :: t                                ! time vector
	real(dp), 		parameter 		:: ti=0._dp, tf=450._dp 			! start and end times
	real(dp) 						:: y1_0,y2_0,y3_0,y4_0 				! initial conditions
	real(dp) 						:: h 								! step
	integer(sp) 					:: i,index_prev						! loop variables
	integer(sp) 					:: istat 							! control variable
    real(dp)                        :: a,b,c                            ! parameters 
    real(dp)                        :: kin_energy,pot_energy,tot_energy ! kinetic, potential and total energies

    a = (1._dp/3._dp) 	! alpha = m2/m1
	b = 0.5_dp		 	! beta  = l2/l1
	c = 0.5_dp			! gamma = g/l1

	y1_0 = 0._dp
	y2_0 = 0._dp
	
    !y3_0 = 0.332_dp
    !y4_0 = 0.845_dp
    
    y3_0 = sqrt(1.125_dp)
	y4_0 = 0._dp

	!subroutine RK4_four_eq(RK4_type, n, a, b, y1_0, y2_0, y3_0,y4_0, y1_RK4, y2_RK4, y3_RK4, y4_RK4, function_type, input_type)
	call RK4_four_eq(1_sp,n,ti,tf,y1_0,y2_0,y3_0,y4_0,y1_RK4,y2_RK4,y3_RK4,y4_RK4,1_sp,1_sp)

    !subroutine energy_dble_pendulum (theta1, theta2, w1,w2, a, b, c, T, U)
    call energy_dble_pendulum (y1_RK4(1), y2_RK4(1), y3_RK4(1),y4_RK4(1), a, b, c, kin_energy, pot_energy)
    tot_energy = (kin_energy + pot_energy)

    t(1) = ti

    open( 10, file = './result.dat', status = 'replace', action = 'write', iostat = istat )
	write(*,*) 'Input/Output file. istat = ', istat
    write(10,*) t(1), y1_RK4(1),y2_RK4(1),y3_RK4(1),y4_RK4(1),tot_energy

    h = abs( tf - ti ) * ( 1._dp / ( real(n,dp) - 1._dp ) )
	do i= 2, n
        index_prev = i-1
		t(i) = t(index_prev) + h

        call energy_dble_pendulum (y1_RK4(i), y2_RK4(i), y3_RK4(i),y4_RK4(i), a, b, c, kin_energy, pot_energy)
        tot_energy = (kin_energy + pot_energy)
        
        write(10,*) t(i), y1_RK4(i),y2_RK4(i),y3_RK4(i),y4_RK4(i),tot_energy
	end do

	close(10)
end program double_pendulum


! OJO! NOTAR QUE LOS RESULTADOS VARIAN MUCHO RESPECTO A LA ELECCIÓN DEL NÚMERO "n"
! La variable tot_energy está en unidades de m1*l1*l1 que son parámetros del primer péndulo