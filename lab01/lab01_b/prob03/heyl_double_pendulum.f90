program heyl_double_pendulum
    use module_presition
    use module_EDO_segundo_orden_flip

    implicit none

    real(dp),       allocatable :: theta1_0(:), theta2_0(:)
    real(dp),       parameter   :: omega1_0=0._dp, omega2_0=0._dp
    real(dp),       parameter   :: ti=0._dp, tf=1000._dp
    integer(sp),    parameter   :: dim_theta1_0=300, dim_theta2_0=150
    real(dp),       parameter   :: theta1_0_min=-3._dp, theta1_0_max=3._dp
    real(dp),       parameter   :: theta2_0_min=0._dp, theta2_0_max=3._dp
    integer(sp)                 :: i,j,istat
    real(dp),       allocatable :: t_flip(:,:)
    integer(sp),    parameter   :: n=50000 ! points number to use in RK4 method
    real(dp)                    :: h_theta1_0,h_theta2_0

    h_theta1_0 = abs(theta1_0_max-theta1_0_min)*(1._dp/(real(dim_theta1_0,dp)-1._dp))
    h_theta2_0 = abs(theta2_0_max-theta2_0_min)*(1._dp/(real(dim_theta2_0,dp)-1._dp))

    ! rellenamos los vectores de condiciones iniciales

    allocate(theta1_0(dim_theta1_0))
    do i=1,dim_theta1_0
        theta1_0(i) = theta1_0_min + h_theta1_0*(real(i,dp)-1._dp)
    end do

    allocate(theta2_0(dim_theta2_0))
    do i=1,dim_theta2_0
        theta2_0(i) = theta2_0_min + h_theta2_0*(real(i,dp)-1._dp)
    end do

    ! Remember
    ! subroutine RK4_four_eq_flip(RK4_type,n,a,b,y1_0,y2_0,y3_0,y4_0,
    !     y1_RK4,y2_RK4,y3_RK4,y4_RK4,x_flip,function_type,input_type)

    ! rellenamos la matriz flip (recordar FORTRAN es colum-major)
    allocate(t_flip(dim_theta2_0,dim_theta1_0))

    open( 10, file = './result_flips.dat', status = 'replace', action = 'write', iostat = istat )
    !write(*,*) 'istat of result_flips.dat = ', istat
    do j=1,dim_theta1_0
        do i=1,dim_theta2_0
            if (i==1_sp .and. j == 1_sp) then
                call RK4_four_eq_flip(n,ti,tf,theta1_0(j),theta2_0(i),omega1_0,omega2_0,t_flip(i,j),.false.)
                write(10,*) theta1_0(j), theta2_0(i), t_flip(i,j)
            else
                call RK4_four_eq_flip(n,ti,tf,theta1_0(j),theta2_0(i),omega1_0,omega2_0,t_flip(i,j),.false.)
                write(10,*) theta1_0(j), theta2_0(i), t_flip(i,j)
            end if
        end do
    end do
    close(10)

! OTRA FORMA DE ESCRIBIRLO, PERO NO ES ÚTIL PARA GNUPLOT
!    open( 10, file = './result_flips.dat', status = 'replace', action = 'write', iostat = istat )
!	 write(*,*) 'istat of result_flips.dat = ', istat
!    do i=1,dim_theta2_0
!        write(10,*) t_flip(i,:)
!    end do
!    close(10)

    deallocate(theta1_0,theta2_0)
    deallocate(t_flip)
    
end program heyl_double_pendulum