!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
! Calculo de la transformada discreta de fourier usando
!  libreria fftw3
!----------------------------------------------------------

program fftw3_b_incise

	use module_presition
	use module_fourier_transform
	
	implicit none
	
	integer(sp) 							:: N 			! discret points number
	integer(sp) 							:: M	 		! nomber of parameters and initial conditions
	real(dp), dimension(:), allocatable 	:: x			! variable vector
	real(dp), dimension(:), allocatable 	:: r			! parameter vector
	real(dp), dimension(:), allocatable		:: omega		! angular frecuency
	real(dp), dimension(:), allocatable		:: sde			! spectral density energy
	real(dp), parameter 					:: x0=0.6_dp	! initial condition
	integer(sp) 							:: i, j, k 		! loop variables
	integer(sp) 							:: istat 		! control variable
	real(dp), dimension(:), allocatable		:: ft_real		! fourier transform real part
	real(dp), dimension(:), allocatable 	:: ft_aimag 	! fourier transform aimag part
	character(len=10)						:: file_id
	character(len=30) 						:: filename

	
	N = 2**9
	M = 5
	
	allocate(r(M))
	r = [1.5_dp, 3.3_dp,3.5_dp,3.55_dp,4._dp]
	
	allocate(x(N), omega(N+1), ft_real(N+1), ft_aimag(N+1), sde(N+1))
	
	20 format (E12.4, E12.4, E12.4, E12.4, E12.4)
	
	x(1) = x0
	do i = 1, M
		do j= 2, N
			x(j) = r(i)*x(j-1)*( 1 - x(j-1) )
		end do
		
		! fft_forware_r2c(N, t_start, t_step, func, ft, sde, omega)
		call fft_forware_r2c(N, 0._dp, 1._dp, x, ft_real,ft_aimag, sde, omega)
		
		! Write the integer into a string:
		write(file_id, '(i0)') i
		! Construct the filename:
		filename = './result_ft_'//trim(adjustl(file_id))//'.dat'
		! Open the file with this name
		open( 10, file = trim(filename), status = 'replace', action = 'write', iostat = istat )
		
		do k = -N/2, N/2
			if ( k < 0 ) then
				write(10,20) omega(-k), ft_real(-k), ft_aimag(-k), sde(-k), r(i)
			else
				write(10,20) omega(N/2+1+k), ft_real(N/2+1+k), ft_aimag(N/2+1+k), sde(N/2+1+k), r(i)
			endif
		enddo
		close(10)
		
	end do
	
	deallocate(r,x,omega,ft_real, ft_aimag,sde)
	
end program fftw3_b_incise
