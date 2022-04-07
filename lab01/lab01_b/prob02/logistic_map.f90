!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
! Calculo de la transformada discreta de fourier usando
!  libreria fftw3
!----------------------------------------------------------

program logistic_map

	use module_presition
	
	implicit none
	
	integer(sp) 							:: N, M, P
	real(dp), dimension(:), allocatable 	:: x, r, x0
	integer(sp) 							:: i, j, k, istat
	
	N = 2**9
	M = 5
	P = 4
	
	
	allocate(r(M))
	r = [1.5_dp, 3.3_dp,3.5_dp,3.55_dp,4._dp]
	allocate(x0(P))
	x0 = [0.06_dp, 0.3_dp,0.6_dp,0.9_dp]
	
	open( 10, file = './result_logistic_map.dat', status = 'replace', action = 'write', iostat = istat )
	20 format (I12, E14.6, E14.6, E14.6)
	
	allocate(x(N))
	
	do i = 1, M
		do k= 1, P
			x(1) = x0(k)
			do j= 1, N
				x(j+1) = r(i)*x(j)*( 1 - x(j) )
				write(10,20) j, x(j), x(1), r(i)
			end do
		end do
	end do
	
	deallocate(x)
	
	close(10)
	
	deallocate(r,x0)
	
end program logistic_map
