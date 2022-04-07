!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
! Calculo de la transformada discreta de fourier usando
!  libreria fftw3
!----------------------------------------------------------

program discret_fourier_transform_fft

	use module_presition
	use module_functions_1D
	
	use, intrinsic :: iso_c_binding 	! Say to Fortran to use C-code
	
	implicit none
	
	include "/usr/include/fftw3.f03" 	! verficar que exista
	
	integer(sp)												:: N 				! number of points for DFT
	type(C_ptr) 											:: plan_rc,plan_cr 	! used plan names of pointer type form C
	integer(sp)	 											:: i, istat 		! loop variable and control variable
	real(dp)												:: t_start, t_end 	! start and end time of temporal window
	real(dp) 												:: t_step, T_samp	! time step and total sampling time
	real(dp) 												:: factor
	real(dp) 												:: omega 						! angular frecuency
	real(dp), parameter 									:: pi = 4._dp*atan(1._dp)
	real(dp), dimension(:), allocatable 					:: t 							! time vector
	real(C_double), dimension(:), allocatable 				:: in_01, out_02
	complex(C_double_complex), dimension(:), allocatable 	:: out_01, in_02
	complex(dp)												:: ft 								! fourier transform
	real(dp) 												:: ft_inverse 						! inverse fourier transform
	real(dp) 												:: sde 								! spectral density energy
	
	write(*,*) 'Input points number for FFT.'
	read(*,*)  N
	write(*,*) 'Input initial sampling time.'
	read(*,*)  t_start
	write(*,*) 'Input final sampling time.'
	read(*,*)  t_end
	
	allocate(in_01(N), out_02(N), out_01(N/2 + 1), in_02(N/2 + 1))
	
	plan_rc = fftw_plan_dft_r2c_1d(N , in_01, out_01, FFTW_MEASURE) 	! create plan r2c
	plan_cr = fftw_plan_dft_c2r_1d(N , in_02, out_02, FFTW_MEASURE)		! create plan c2r
	
	! Function to transform
	
	T_samp = abs(t_end-t_start)
	t_step = T_samp*(1._dp/real(N,dp))

	open( 11, file = './result_function.dat', status = 'replace', action = 'write', iostat = istat )
	21 format (E12.4, E12.4)
	if (istat == 0_sp) then
		allocate(t(N))
		do i = 1, N
			t(i) = t_start + t_step*(real(i,dp)-1._dp)
			
			!in_01(i) = exp(-t(i)**2)*t(i) 													! function_01
			!in_01(i) = exp(-pi*(t(i)**2))*cos(20._dp*pi*t(i)) 								! function_02
			
			! f1D_trig(x,a1,a2,p1,p2,c,function_type) => cos(a1*(x**p1))+sin(a2*(x**p2))+c
			in_01(i) = f1D_trig(t(i),(20._dp*pi), (pi*0.5_dp), 1._dp, 1._dp, 0._dp,1_sp) 	! function_03
			out_02(i) = in_01(i)
			write(11,21) t(i), in_01(i)
		end do
	else
		write(*,*) 'Error opening file. iostat != 0'
	end if
	
	!----------------------------------------------------------
	! COMMENTS
	!----------------------------------------------------------
	! One-Dimensional DFTs of Real Data
	!  "r2c": the real input to complex-Hermitian output transform (always FFTW_FORWARD)
	!  "c2r": complex-Hermitian input to real output transform (always FFTW_BACKWARD)
	!
	!  Fourier transform
	!    fftw_plan fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags)
	!  Inverse fourier transform
	!    fftw_plan fftw_plan_dft_c2r_1d(int n, fftw_complex *in, double *out, unsigned flags)
	!
	!  plan_rc = fftw_plan_dft_r2c_1d(N , in , out, flags) ! create plan r2c
	!  plan_cr = fftw_plan_dft_c2r_1d(N , in , out, flags) ! create plan c2r
	!  flags: [FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT]
	!----------------------------------------------------------
	
	! Fourier transform calculus
	call fftw_execute_dft_r2c(plan_rc, in_01 , out_01) 	! run forward transform

	open( 10, file = './result_transforms.dat', status = 'replace', action = 'write', iostat = istat )
	20 format (E12.4, E12.4, E12.4, E12.4)
	if (istat == 0_sp) then
		factor = 1._dp/T_samp
		do i = -N/2, N/2
			omega = real(i,dp)*2._dp*pi*factor
			if ( i < 0 ) then
				ft = conjg(out_01(-i+1))/real(N,dp)
				in_02(-i) = ft
				sde = ft*conjg(ft) ! or sde = abs(ft)**2
				write(10,20) omega, real(ft), aimag(ft), sde
			else
				ft = out_01(i + 1)/real(N,dp)
				in_02(i) = ft
				sde = ft*conjg(ft) ! or sde = abs(ft)**2
				write(10,20) omega, real(ft), aimag(ft), sde
			endif
		enddo
	else
		write(*,*) 'Error opening file. iostat != 0'
	end if
	close(10)
	
	! Inverse fourier transform calculus
	call fftw_execute_dft_c2r(plan_rc,in_02,out_02) ! run backward transform
	
	open( 12, file = './result_inverse_transforms.dat', status = 'replace', action = 'write', iostat = istat )
	if (istat == 0_sp) then
		do i = 1, N
			ft_inverse = out_02(i)
			write(12,21) t(i), ft_inverse
		end do
	else
		write(*,*) 'Error opening file. iostat != 0'
	end if
	close(12)
	
	deallocate(t)
	
	call fftw_destroy_plan(plan_rc)
	call fftw_destroy_plan(plan_cr)
	
end program discret_fourier_transform_fft


!----------------------------------------------------------
! REFERENCES
!----------------------------------------------------------
! https://stackoverflow.com/questions/44626196/fortran-and-fftw3
! https://www.fftw.org/fftw3.pdf
! https://www-uxsup.csx.cam.ac.uk/courses/moved.Fortran/paper_07.pdf
! https://www.gnu.org/software/libc/manual/html_node/Absolute-Value.html
! https://www.fftw.org/fftw3_doc/Fortran-Examples.html
!----------------------------------------------------------

!----------------------------------------------------------
! COMPILATION RUL
!----------------------------------------------------------
! gfortran -o test.x  test.f90  -lfftw3
!----------------------------------------------------------
