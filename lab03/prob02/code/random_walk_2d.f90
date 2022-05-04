! Problema 02
program random_walk_2d
    use module_precision

    implicit none
    integer(sp)            :: seed     ! semilla
    integer(sp), parameter :: n=1E6_sp ! numero total de random walks
    real(dp)               :: x,y      ! pasos
    integer(sp), parameter :: n_step=1E+03_sp ! numero de pasos de c/ random walk
    real(dp)               :: cuad01,cuad02,cuad03,cuad04 ! contadores en c/ cuadrante
    real(dp)               :: cuad01_tot,cuad02_tot,cuad03_tot,cuad04_tot ! contadores en c/ cuadrante
    integer(sp)            :: i
    

    cuad01_tot=0._dp;cuad02_tot=0._dp;cuad03_tot=0._dp;cuad04_tot=0._dp
    do i=1,n
        call walk_2d(0_sp,seed,1_sp,n_step,0._dp,0._dp,x,y,cuad01,cuad02,cuad03,cuad04)
        cuad01_tot=cuad01_tot+cuad01
        cuad02_tot=cuad02_tot+cuad02
        cuad03_tot=cuad03_tot+cuad03
        cuad04_tot=cuad04_tot+cuad04
    end do

    cuad01=cuad01_tot*(1._dp/n)
    cuad02=cuad02_tot*(1._dp/n)
    cuad03=cuad03_tot*(1._dp/n)
    cuad04=cuad04_tot*(1._dp/n)
    
    write(*,'(A10,E10.2)') 'N/4=', real(n_step)*0.25_dp
    write(*,'(A10,E10.2)') 'cuad01=', cuad01
    write(*,'(A10,E10.2)') 'cuad02=', cuad02
    write(*,'(A10,E10.2)') 'cuad03=', cuad03
    write(*,'(A10,E10.2)') 'cuad04=', cuad04
    write(*,'(A10,I10)') 'N=', n_step
    write(*,'(A10,E10.2)') 'Naprox=', cuad01+cuad02+cuad03+cuad04


end program random_walk_2d

! subrutina para realizar caminata aleatoria de n_step pasos específicos
! partiendo de cierto origen especificoe, elegir tipo pseudo generador random,
! devolver, semilla y contar cuantos pasos cayeron en determinado cuadrante
subroutine walk_2d(switch,seed,rnd_type,n_step,x0,y0,x,y,&
    count_cuad_01,count_cuad_02,count_cuad_03,count_cuad_04)
    use module_precision;use module_random_generator
    use module_mzran;use module_mt19937

    implicit none
    integer(sp), intent(in)    :: switch ! prender o apagar escritura de datos
    integer(sp), intent(inout) :: seed
    integer(sp), intent(in)    :: rnd_type,n_step
    real(dp), intent(in)       :: x0,y0
    real(dp), intent(inout)    :: x,y
    real(dp), intent(inout)    :: count_cuad_01,count_cuad_02,&
                                  count_cuad_03,count_cuad_04

    real(dp),    parameter :: px=0.5_dp,py=0.5_dp ! probabilidades de pasos
    real(dp),    parameter :: step=1._dp          ! longitud de paso (fija)
    integer(sp)            :: seed_val(8),i,istat
    real(dp)               :: nrand

    ! generamos la semilla
    call date_and_time(values=seed_val)
    seed=seed_val(8)*seed_val(7)*seed_val(6)+seed_val(5)

    if (switch==1_sp) then
        open(10,file='../results/result.dat',status='replace',action='write',iostat=istat)
        20 format(E11.4,x,E11.4)
        write(*,*) 'istat=', istat
    end if
    ! posición inicial
    x=x0;y=y0; if (switch==1_sp) write(10,20) x,y
    count_cuad_01=0._dp;count_cuad_02=0._dp;count_cuad_03=0._dp;count_cuad_04=0._dp

    select case(rnd_type)
        case(1);nrand=ran0(seed)                       ! ran0 random generator
        case(2);nrand=ran2(seed)                       ! ran2 random generator
        case(3);nrand=rmzran()                         ! mzran random generator
        case(4);call sgrnd(seed);nrand=real(grnd(),dp) ! mt19937 random generator
        case default; write(*,*) 'Invalid random generator type'
    end select

    do i=2,n_step
        nrand=ran0(seed)
        cond1:   if (0._dp<nrand.and.nrand<0.25_dp)  then; x=x+step; exit cond1
            else if (0.25_dp<nrand.and.nrand<0.5_dp) then; y=y+step; exit cond1
            else if (0.5_dp<nrand.and.nrand<0.75_dp) then; x=x-step; exit cond1
            else if (0.75_dp<nrand.and.nrand<1._dp)  then; y=y-step; exit cond1
        end if cond1
        if (switch==1_sp) write(10,20) x,y
        ! Determinamos el cuadrante de la partícula
        cond2:   if (x>0._dp.and.y>0._dp) then; count_cuad_01=count_cuad_01+1._dp; exit cond2
            else if (x<0._dp.and.y>0._dp) then; count_cuad_02=count_cuad_02+1._dp; exit cond2
            else if (x<0._dp.and.y<0._dp) then; count_cuad_03=count_cuad_03+1._dp; exit cond2
            else if (x>0._dp.and.y<0._dp) then; count_cuad_04=count_cuad_04+1._dp; exit cond2
            else if (x==0._dp.and.y==0._dp) then     ! origen de coordenadas
            count_cuad_01=count_cuad_01+0.25_dp;count_cuad_02=count_cuad_02+0.25_dp
            count_cuad_03=count_cuad_03+0.25_dp;count_cuad_04=count_cuad_04+0.25_dp; exit cond2
            else if (x>0._dp.and.y==0._dp) then ! semi-eje x positivo
            count_cuad_01=count_cuad_01+0.5_dp;count_cuad_04=count_cuad_04+0.5_dp; exit cond2
            else if (x<0._dp.and.y==0._dp) then ! semi-eje x negativo
            count_cuad_02=count_cuad_02+0.5_dp;count_cuad_03=count_cuad_03+0.5_dp; exit cond2
            else if (x==0._dp.and.y>0._dp) then ! semi-eje y positivo
            count_cuad_01=count_cuad_01+0.5_dp;count_cuad_02=count_cuad_02+0.5_dp; exit cond2
            else if (x==0._dp.and.y<0._dp) then ! semi-eje y negativo
            count_cuad_03=count_cuad_03+0.5_dp;count_cuad_04=count_cuad_04+0.5_dp; exit cond2
        end if cond2
    end do
    if (switch==1_sp) close(10)
end subroutine walk_2d