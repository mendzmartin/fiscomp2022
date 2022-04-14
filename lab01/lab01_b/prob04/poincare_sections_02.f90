program poincare_sections_02
    use module_presition
    use module_EDO_segundo_orden_poincare
    use module_pullen_edmonds

    implicit none
    real(dp), parameter     :: q2_0=0._dp               ! coordenada generaliza 2 (cond. inicial)
    integer(sp), parameter  :: n=1500000_sp             ! puntos para integrar con RK4
    integer(sp), parameter  :: nloop = 18_sp            ! Número de diferentes condiciones iniciales
    integer(sp), parameter  :: np = 20000_sp            ! Número de puntos intersección a grabar
    real(dp), parameter     :: ti=0._dp,tf=15000._dp    ! tiempos de integración
    integer(sp)             :: i,j,k,l                  ! variables de loop
    real(dp)                :: HE,m0,m1
    real(dp)                :: p2_0                     ! momento generalizado 2 (cond. inicial)
    real(dp), allocatable   :: E(:)
    real(dp), allocatable   :: p1_0(:),q1_0(:)          ! coord. y mom. generalizados 1 (cond. inicial)
    real(dp), allocatable   :: p1(:),q1(:)              ! coord. y mom. generalizados 1 (sol. RK4)
    character(len=10)       :: file_id
	character(len=30)       :: filename

    ! Vector de energías
    E = [5._dp,20._dp,100._dp]

    allocate(q1_0(nloop),p1_0(nloop))
    allocate(q1(np),p1(np))
    do1: do i = 1_sp,3_sp
        write(file_id, '(i0)') i
        filename = './resultados_E'//trim(adjustl(file_id))//'.dat'
        open( 10, file = trim(filename), status = 'replace', action = 'write')
        ! Calculo las condiciones iniciales
        m0 = sqrt(2._dp*E(i))
        p1_0 = -m0+2._dp*m0*(/(real(j,dp), j = 0,nloop)/)*(1._dp/real(nloop,dp))
        q1_0 = -m0+2._dp*m0*(/(real(j,dp), j = 0,nloop)/)*(1._dp/real(nloop,dp))
        do2: do j = 1,nloop,2_sp
            do3: do k = 1,nloop,2
                m1 = 2._dp*E(i)-q1_0(k)*q1_0(k)-p1_0(j)*p1_0(j)
                ! Restricción de p1 y q1
                if (m1 > 0._dp) then
                    ! Defino p2_0
                    p2_0 = sqrt(m1)
                    !energy_pullen_edmons(q1,q2,p1,p2,alpha,m,omega,E,function_type)
                    call energy_pullen_edmons(q1_0(k),q2_0,p1_0(j),p2_0,0.05_dp,1._dp,1._dp,HE,1_sp)
                    ! Integro el sistema en el tiempo con la subrutina rk4_PEH_Poincare
                    call RK4_four_eq_poincare(n,np,ti,tf,q1_0(k),q2_0,p1_0(j),p2_0,p1,q1)
                    do4: do l = 1,np
                        ! Escribo valores
                        write(10,*) p1(l),q1(l)
                    enddo do4
                endif
            enddo do3
        enddo do2
        close(10)
    enddo do1

end program poincare_sections_02