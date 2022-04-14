program pullen_edmonds
    
    use module_presition
    use module_EDO_segundo_orden

    implicit none
    integer(sp), parameter :: n_opt=4194641 ! aprox 2^22 (n_opt = (8/eps_machine)^(2/5))
    real(dp), parameter :: ti=0._dp,tf=200._dp
    real(dp), parameter :: q1_0=2._dp,q2_0=0._dp,p1_0=0._dp
    real(dp) :: p2_0
    real(dp), allocatable :: E(:),q1_RK4(:),q2_RK4(:),p1_RK4(:),p2_RK4(:)
    character(len=10)						:: file_id
	character(len=30) 						:: filename
    integer(sp) :: i,j,istat

    allocate(E(3))
    E=[5._dp,20._dp,100._dp]

    allocate(q1_RK4(n_opt),q2_RK4(n_opt),p1_RK4(n_opt),p2_RK4(n_opt))
    do i=1,size(E)
        p2_0=sqrt(2._dp*E(i)-4._sp)

        ! RK4_four_eq(RK4_type,n,a1,b1,y1_0,y2_0,y3_0,y4_0,y1_RK4,y2_RK4,y3_RK4,y4_RK4,function_type,input_type)
        call RK4_four_eq(1_sp,n_opt,ti,tf,q1_0,q2_0,p1_0,p2_0,q1_RK4,q2_RK4,p1_RK4,p2_RK4,1_sp,2_sp)

		write(file_id, '(i0)') i
		filename = './results_E'//trim(adjustl(file_id))//'.dat'
        open( 10, file = trim(filename), status = 'replace', action = 'write', iostat = istat )

        do j=1,n_opt
            write(10,*) q1_RK4(j),q2_RK4(j),p1_RK4(j),p2_RK4(j)
        end do

        close(10)

    end do

    deallocate(E)
    deallocate(q1_RK4,q2_RK4,p1_RK4,p2_RK4)

end program pullen_edmonds