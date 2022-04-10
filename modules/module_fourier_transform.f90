!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
! Calculo de la transformada discreta de fourier usando
!  libreria fftw3
!----------------------------------------------------------

module module_fourier_transform
	use module_presition
	
	implicit none

	contains
	
	subroutine fft_forware_r2c(N, t_start, t_step, func, ft_real, ft_aimag, sde, omega)
		use module_presition
		use, intrinsic :: iso_c_binding 	! Say to Fortran to use C-code
		
		implicit none
		
		include "/usr/include/fftw3.f03" 	! verficar que exista
		
		integer(sp), 					intent(in)				:: N 				! number of points for DFT
		real(dp), 						intent(in)				:: t_start, t_step 	! start and end time of temporal window
		real(dp), dimension(N+1), 		intent(out)				:: ft_real 				! fourier transform
		real(dp), dimension(N+1), 		intent(out)				:: ft_aimag
		real(dp), dimension(N+1), 		intent(out) 			:: sde 				! spectral density energy
		real(C_double), dimension(N), 	intent(in)				:: func
		real(dp), dimension(N+1), 		intent(out)				:: omega 					! angular frecuency
		
		
		type(C_ptr) 											:: plan_rc		 	! used plan names of pointer type form C
		integer(sp)	 											:: i		! loop variable and control variable
		real(dp) 												:: T_samp	! time step and total sampling time
		real(dp) 												:: factor
		real(dp), parameter 									:: pi=4._dp*atan(1._dp)
		!real(dp), dimension(:), allocatable 					:: t 						! time vector
		real(C_double), dimension(:), allocatable				:: in_01
		complex(C_double_complex), dimension(:), allocatable 	:: out_01
		complex(dp) 											:: ft
		
		allocate(in_01(N), out_01(N/2 + 1))
		plan_rc = fftw_plan_dft_r2c_1d(N , in_01, out_01, FFTW_MEASURE) 	! create plan r2c
		
		! Function to transform
		T_samp = t_step*real(N,dp)

		!allocate(t(N))
		do i = 1, N
			!t(i) = t_start + t_step*(real(i,dp)-1._dp)
			in_01(i) = func(i) 	! function_03
		end do
		
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
		
		factor = 1._dp/T_samp
		do i = -N/2, N/2
			if ( i < 0 ) then
				omega(-i) = real(i,dp)*2._dp*pi*factor
				ft = conjg(out_01(-i+1))/real(N,dp)
				ft_real(-i) = real(ft)
				ft_aimag(-i) = aimag(ft)
				sde(-i) = ft*conjg(ft) ! or sde = abs(ft)**2
			else
				omega(N/2+1+i) = real(i,dp)*2._dp*pi*factor
				ft = out_01(i + 1)/real(N,dp)
				ft_real(N/2+1+i) = real(ft)
				ft_aimag(N/2+1+i) = aimag(ft)
				sde(N/2+1+i) = ft*conjg(ft) ! or sde = abs(ft)**2
			endif
		enddo
		
		!deallocate(t)
		call fftw_destroy_plan(plan_rc)
		
	end subroutine fft_forware_r2c
end module module_fourier_transform


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
