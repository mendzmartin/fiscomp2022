! Problema 02
program random_walk_2d
    use module_precision

    implicit none
    integer(sp)            :: seed,seed_val(8)                            ! semilla
    integer(sp), parameter :: n_step_total=1000000_sp                     ! numero total de random walks
    integer(sp), parameter :: switch=3_sp                                 ! cambiar si se quieren escribir los datos
    integer(sp), parameter :: n_step=1000_sp                              ! numero de pasos de c/ random walk
    real(dp)               :: x,y,x_old,y_old,dcm,dcm_tot                 ! pasos y desplazamiento cuadrático medio
    real(dp)               :: cuad01,cuad02,cuad03,cuad04                 ! contadores en c/ cuadrante
    real(dp)               :: cuad01_tot,cuad02_tot,cuad03_tot,cuad04_tot ! contadores en c/ cuadrante
    integer(sp)            :: i,j,istat
    integer(sp)            :: rnd_type                                    ! tipo de random generator
    real(dp)               :: suma                                        ! variable de control

    open(10,file='../results/result.dat',status='replace',action='write',iostat=istat)
    select case(switch)
        case(1) ! mapa random walk
            open(10,file='../results/result_01.dat',status='replace',action='write',iostat=istat)
            21 format(A12,x,A12); write(10,21) 'x-coord','y-coord'
        case(2) ! inciso a
            open(10,file='../results/result_02.dat',status='replace',action='write',iostat=istat)
            22 format(A12,x,A12); write(10,22) 'n_step','dcm'
            23 format(I12,x,E12.4)
        case(3) ! inciso b
            open(10,file='../results/result_03.dat',status='replace',action='write',iostat=istat)
            24 format(5(A12,x),A12); write(10,24) 'j','N/4','cuad01','cuad02','cuad03','cuad04'
            25 format(I12,x,5(E12.4,x),E12.4)
    end select
    if (istat /= 0_sp) write(*,*) 'istat_error=', istat

    rnd_type=4_sp ! elegir random generator
    do j=100_sp,n_step,100_sp
        write(*,*) j
        ! generamos la semilla
        call date_and_time(values=seed_val)
        seed=seed_val(8)*seed_val(7)*seed_val(6)+seed_val(5)
        x_old=0._dp;y_old=0._dp;dcm_tot=0._dp
        cuad01_tot=0._dp;cuad02_tot=0._dp;cuad03_tot=0._dp;cuad04_tot=0._dp
        do i=1,n_step_total
            call walk_2d(switch,10,rnd_type,seed,j,x_old,y_old,x,y,dcm,cuad01,cuad02,cuad03,cuad04)
            cuad01_tot=cuad01_tot+cuad01;cuad02_tot=cuad02_tot+cuad02
            cuad03_tot=cuad03_tot+cuad03;cuad04_tot=cuad04_tot+cuad04
            dcm_tot=dcm_tot+dcm
            if (switch==1_sp) then;x_old=x;y_old=y; end if
        end do

        if (switch==2_sp) write(10,23) j,dcm_tot*(1._dp/n_step_total)

        cuad01=cuad01_tot;cuad02=cuad02_tot;cuad03=cuad03_tot;cuad04=cuad04_tot
        suma=(cuad01+cuad02+cuad03+cuad04)*0.25_dp

        if (switch==3_sp) write(10,25) j,real(n_step_total)*0.25_dp,cuad01,cuad02,cuad03,cuad04,suma
    end do
    close(10)
end program random_walk_2d

! subrutina para realizar caminata aleatoria de n_step pasos específicos
! partiendo de cierto origen especificoe, elegir tipo pseudo generador random,
! devolver, semilla y contar cuantos pasos cayeron en determinado cuadrante
subroutine walk_2d(switch,n_file,rnd_type,seed,n_step,x0,y0,x,y,dcm,&
    count_cuad_01,count_cuad_02,count_cuad_03,count_cuad_04)
    use module_precision;use module_random_generator
    use module_mzran;use module_mt19937

    implicit none
    integer(sp), intent(in)    :: switch,n_file ! prender o apagar escritura de datos
    integer(sp), intent(in)    :: n_step    ! numero de pasos totales
    real(dp),    intent(in)    :: x0,y0     ! cordenadas iniciales
    integer(sp), intent(inout) :: rnd_type,seed
    real(dp),    intent(out)   :: x,y      ! coordenadas
    real(dp),    intent(out)   :: count_cuad_01,count_cuad_02,&
                                  count_cuad_03,count_cuad_04
    real(dp),   intent(out)    :: dcm       ! desplazamiento cuadrático medio

    !real(dp),    parameter :: px=0.5_dp,py=0.5_dp ! probabilidades de pasos
    real(dp),    parameter :: step=1._dp           ! longitud de paso (fija)
    integer(sp)            :: i
    real(dp)              :: nrand                ! numero pseudo-aleatorio

    if (switch==1_sp) then;20 format(E12.4,x,E12.4);write(n_file,20) x,y;end if
    ! posición inicial
    x=x0;y=y0;count_cuad_01=0._dp;count_cuad_02=0._dp;count_cuad_03=0._dp;count_cuad_04=0._dp
    !if (rnd_type==4_sp) call sgrnd(seed)
    do i=1,n_step
        select case(rnd_type)
                case(1);nrand=ran0(seed)                       ! ran0 random generator
                case(2);nrand=ran2(seed)                       ! ran2 random generator
                case(3);nrand=rmzran()                         ! mzran random generator
                case(4);nrand=real(grnd(),dp) ! mt19937 random generator
                case default; write(*,*) 'Invalid random generator type'
        end select
        cond1:   if (0._dp<=nrand.and.nrand<0.25_dp)  then; x=x+step; exit cond1
            else if (0.25_dp<=nrand.and.nrand<0.5_dp) then; y=y+step; exit cond1
            else if (0.5_dp<=nrand.and.nrand<0.75_dp) then; x=x-step; exit cond1
            else if (0.75_dp<=nrand.and.nrand<1._dp)  then; y=y-step; exit cond1
        end if cond1
        if (switch==1_sp) write(n_file,20) x,y
    end do
    ! Determinamos el cuadrante de la partícula al final de la caminata
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
    if (switch==2_sp) dcm=(x-x0)*(x-x0)+(y-y0)*(y-y0)
end subroutine walk_2d