program poincare_sections
    use module_presition
    use module_EDO_segundo_orden
    use module_pullen_edmonds

    implicit none
    integer(sp), parameter :: n_opt=20000  !4194641 ! aprox 2^22 (n_opt = (8/eps_machine)^(2/5))
    integer(sp), parameter :: n_ic=100
    real(dp), parameter :: ti=0._dp,tf=15000._dp
    real(dp) :: q1_0,q2_0,p1_0
    real(dp), parameter     :: p2_0 = 10_dp
    real(dp), allocatable :: E(:),q1_RK4(:),q2_RK4(:),p1_RK4(:),p2_RK4(:)
    integer(sp) :: i,j,k,l,istat,count
    real(dp) :: hq1,hp1
    real(dp), parameter :: q1_min=-20._dp, q1_max=20._dp
    real(dp), parameter :: p1_min=-20._dp, p1_max=20._dp
    logical :: control

    allocate(E(3))
    E=[5._dp,20._dp,100._dp]

    allocate(q1_RK4(n_opt),q2_RK4(n_opt),p1_RK4(n_opt),p2_RK4(n_opt))
    open( 10, file = 'poincare.dat', status = 'replace', action = 'write', iostat = istat )
    write(*,*) istat
    do i=1,size(E)
        do k=1,n_ic
            hq1 = abs(q1_max - q1_min)*(1._sp/(real(n_ic,dp)-1._dp))
            q1_0=q1_min + hq1*(real(k,dp)-1._dp)
            do l=1,n_ic
                hp1 = abs(p1_max - p1_min)*(1._sp/(real(n_ic,dp)-1._dp))
                p1_0=p1_min + hp1*(real(l,dp)-1._dp)

                ! q2_pullen_edmons(q1,q2,p1,p2,alpha,m,omega,E)
                call q2_pullen_edmons(q1_0,q2_0,p1_0,p2_0,0.05_dp,1._dp,1._dp,E(i),control)

                if (control .eqv. .true. ) then
                    ! RK4_four_eq(RK4_type,n,a1,b1,y1_0,y2_0,y3_0,y4_0,y1_RK4,y2_RK4,y3_RK4,y4_RK4,function_type,input_type)
                    call RK4_four_eq(1_sp,n_opt,ti,tf,q1_0,q2_0,p1_0,p2_0,q1_RK4,q2_RK4,p1_RK4,p2_RK4,1_sp,2_sp)

                    count = 0_sp
                    do j=2,n_opt
                        if (q2_RK4(j)*q2_RK4(j-1_sp) <= 0._dp)  then
                            count=count+1_sp
                            if (count < 50_sp) then
                                q1_RK4(j)=0.5_dp*(q1_RK4(j)+q1_RK4(j-1_sp))
                                p1_RK4(j)=0.5_dp*(p1_RK4(j)+p1_RK4(j-1_sp))
                                write(*,*) q1_RK4(j),p1_RK4(j)
                            end if
                        else
                            write(*,*) q1_RK4(j),p1_RK4(j)
                        end if
                    end do
                else
                    write(*,*) 'error'
                end if
            end do
        end do
    end do
    close(10)
    deallocate(q1_RK4,q2_RK4,p1_RK4,p2_RK4)
    deallocate(E)

end program poincare_sections
