!----------------------------------------------------------
! PURPOSE
!----------------------------------------------------------
! Points and weights for Gaussian quadrature
! Rescale rescales the gauss-legendre grid points and weights
! npts	:= number of points
! job	:=	0 rescalling uniformly between (a,b)
! 			1 for integral (0,b) with 50% inside (a,b+2a)
! 			2 for integral (a,inf) with 50% inside (a,b+2a)
! x, w 	:= output grid points and weights
!----------------------------------------------------------

module module_gauss

	use module_presition, only: pr => dp, sp
	
	contains
	
	subroutine gauss(npts, job, a, b, x, w)
		implicit none
		
		integer(sp), intent(in) :: npts, job
		real(pr), intent(in) 	:: a, b
		real(pr), intent(out) 	:: x(0:npts),w(0:npts)
		real(pr) 				:: t, t1, pp, p1, p2, p3, eps
		real(pr) 				:: pi
		real(pr) 				:: xi
		integer(sp)				:: m, i, j

		pi = acos(-1._pr)
		m = 0; i = 0; j = 0
		t = 0._pr; t1 = 0._pr; pp = 0._pr; 
		p1 = 0._pr; p2 = 0._pr; p3 = 0._pr
		eps = 3.e-14_pr 					! Accuracy: ********ADJUST THIS*********!
		m = (npts + 1)/2
!		do i = 1, (m + 1)
		do i = 1, (m + 1)
!			t = cos(pi*(dble(i) - 0.25_pr)/(dble(npts) + 0.5_pr) )
			t = cos(pi*(dble(i) - 0.25_pr)/(dble(npts+1) + 0.5_pr))
			t1 = 1._pr
			do while( (abs(t - t1) ) >= eps)
				p1 = 1._pr ; p2 = 0._pr
				
				do j=1,npts + 1
					p3 = p2; p2 = p1
! 					p1 = ( (2._pr*dble(j) - 1._pr)*t*p2 - (dble(j) - 1._pr)*p3)/(dble(j) )
					p1 = ( (2._pr*dble(j) - 1._pr)*t*p2 - dble(j-1)*p3)/(dble(j) )
				enddo
				
!				pp = dble(npts)*(t*p1 - p2)/(t*t - 1._pr)
				pp = dble(npts+1)*(t*p1 - p2)/(t*t - 1._pr)
				t1 = t
				t = t1 - p1/pp
			enddo
			
			x(i - 1) = - t; 
! 			x(npts - i) = t
			x(npts + 1 - i) = t
			w(i - 1) = 2._pr/( (1._pr - t*t)*pp*pp)
! 			w(npts - i) = w(i - 1)
			w(npts + 1 - i) = w(i - 1)
!			print *," x(i - 1)", x(i - 1) , " w " , w(npts - i),i

		enddo
		
		! Rescale the grid points
		
		! Scale to (a,b) uniformly
		if (job == 0) then
				do i=0, npts
					x(i) = x(i)*(b - a)/2._pr + (b + a)/2._pr
					w(i) = w(i)*(b - a)/2._pr
				enddo
			
			! Scale to (0,b) with 50% pints inside (0,ab/(a+b))
		else if (job == 1) then
			do i=0, npts
				xi = x(i)
				x(i) = a*b*(1._pr + xi) / (b + a - (b - a)*xi)
				w(i) = w(i)*2._pr*a*b*b/( (b + a - (b - a)*xi)*(b + a - (b - a)*xi) )
			enddo
		
		! Scale to (a,inf) with 50% points inside (a,b+2a)
		else if (job == 2) then
			do i=0, npts
				xi = x(i)
				x(i) = (b*xi + b + a + a) / (1._pr - xi)
				w(i) = w(i)*2.*(a + b)/( (1._pr - xi)*(1. - xi) )
			enddo
		
		else
			write(*,*) "Wrong value of job"
		endif

		return
		
	end subroutine gauss
	
end module module_gauss
